#!/usr/bin/env python
"""Drop-in wrapper around `cellbender remove-background` that fixes the
checkpoint-save crash.

Problem
-------
CellBender 0.3.x checkpoints the model at the final epoch via
`torch.save(model_obj, ...)` (remove_background/checkpoint.py:save_checkpoint).
pyro registers the model's parameters with `pyro.module(..., update_module_params=True)`,
which attaches a `.unconstrained` attribute to every `nn.Parameter` for its
constraint bookkeeping -- and that attribute is a `weakref.ReferenceType`.
`torch.save` (pickle) cannot serialize a weakref, so the save aborts with
    TypeError: cannot pickle 'weakref.ReferenceType' object
No `ckpt.tar.gz` is written, and the subsequent posterior step then dies with
    AssertionError: Checkpoint file ckpt.tar.gz does not exist ...

This is independent of torch and pyro versions (reproduced on torch 2.3.1 / 2.5.1 /
2.13 and pyro 1.8.6 / 1.9.1), because the weakref comes from pyro's constraint
mechanism, not from any version-specific scheduler/optimizer behaviour. Upstream
issues broadinstitute/CellBender #395 / #386 / #296 are the same bug, unresolved.

Fix
---
Monkeypatch `save_checkpoint` so that, immediately before it serializes the model,
we strip the `.unconstrained` weakref off every reachable Parameter, let the save
proceed, then restore the attributes. Note the CompositeEncoder is a plain `dict`
(not a registered nn.Module), so `model.parameters()` does NOT reach the encoder's
params -- we recurse the object graph explicitly to catch them.

Usage: identical to the cellbender CLI, e.g.
    python run_cellbender_patched.py remove-background --input in.h5 --output out.h5 ...
Must run inside the cellbender conda env (so `import cellbender` resolves).
"""
import sys

import torch
import cellbender.remove_background.checkpoint as _ck


def _collect_params(root, seen=None, out=None, depth=0):
    """Recursively gather every nn.Parameter reachable from `root`, walking into
    nn.Modules, dicts, lists/tuples/sets, and plain objects' __dict__ -- so the
    CompositeEncoder dict's sub-encoders are covered too."""
    if seen is None:
        seen, out = set(), []
    if id(root) in seen or depth > 10:
        return out
    seen.add(id(root))
    if isinstance(root, torch.nn.Parameter):
        out.append(root)
    if isinstance(root, torch.nn.Module):
        children = (list(root._parameters.values()) + list(root._modules.values())
                    + list(root._buffers.values()) + list(vars(root).values()))
    elif isinstance(root, dict):
        children = list(root.values())
    elif isinstance(root, (list, tuple, set)):
        children = list(root)
    elif hasattr(root, '__dict__'):
        children = list(vars(root).values())
    else:
        children = []
    for c in children:
        if c is not None:
            _collect_params(c, seen, out, depth + 1)
    return out


_orig_save_checkpoint = _ck.save_checkpoint


def _patched_save_checkpoint(filebase, model_obj, scheduler, *args, **kwargs):
    stripped = []
    for p in _collect_params(model_obj):
        if hasattr(p, 'unconstrained'):
            stripped.append((p, p.unconstrained))
            try:
                object.__delattr__(p, 'unconstrained')
            except Exception:
                try:
                    delattr(p, 'unconstrained')
                except Exception:
                    pass
    try:
        return _orig_save_checkpoint(filebase, model_obj, scheduler, *args, **kwargs)
    finally:
        # restore pyro's bookkeeping so in-memory training/posterior is unaffected
        for p, val in stripped:
            try:
                p.unconstrained = val
            except Exception:
                pass


_ck.save_checkpoint = _patched_save_checkpoint

if __name__ == '__main__':
    from cellbender.base_cli import main
    sys.exit(main())
