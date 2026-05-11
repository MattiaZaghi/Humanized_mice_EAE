library(ShinyCell2)
library(Seurat)
#in R

hsmm_ozge <- readRDS("/date/gcb/gcb_MZ/Humanized_mice_OD/Results/objects/humanized_mouse_spleen_whole_doublets.rds")

# ShinyCell2 requires a single joined RNA layer as default assay
DefaultAssay(hsmm_ozge) <- "RNA"
hsmm_ozge[["RNA"]] <- JoinLayers(hsmm_ozge[["RNA"]])

scConf <- createConfig(hsmm_ozge)

makeShinyFiles(hsmm_ozge, scConf, shiny.prefix="Spleen", shiny.dir="/date/gcb/gcb_MZ/Humanized_mice_OD/Results/shinyapps/Spleen")
makeShinyCodes(shiny.title = "Spleen", shiny.prefix="Spleen",
               shiny.dir="/date/gcb/gcb_MZ/Humanized_mice_OD/Results/shinyapps/Spleen")


# once that runs if you want to host it in our shinny account
#these are the credentials
#remenber to remove this part if you would share this code

rsconnect::setAccountInfo(name='castelobranco',
                          token='158254637BEBE253C782F8B4D7F38F40',
                          secret='Pg+JjWA9Gt+Lf5FJlO6QiyxY9GRXE1JHBy+CpVpM')

user <- "castelobranco" 

rsconnect::deployApp("/date/gcb/gcb_MZ/Humanized_mice_OD/Results/shinyapps/Spleen",    account=user,      server="shinyapps.io")

# when this finishes a window open and you just need to click yes or "try again"

#if you would create any other type of shiny app you can always launch them with the secret credentials. Be careful to not share these ones
