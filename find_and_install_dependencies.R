library(packrat)


install_missing_pkgs <- function(path, do_install = TRUE){
  pkgs <- packrat:::appDependencies(path)
  installed_pkgs <- installed.packages()
  cran_pkgs <- available.packages()[, "Package"]

  if(!"BiocManager" %in% installed_pkgs) {
    install.packages("BiocManger")
  }
  bioc_pkgs <- setdiff(BiocManager::available(include_installed = FALSE),
                       cran_pkgs)

  new_packages <- pkgs[!(pkgs %in% installed_pkgs[,"Package"])]

  needed_cran_pkgs <- new_packages[new_packages %in% cran_pkgs]
  needed_bioc_pkgs <- new_packages[new_packages %in% bioc_pkgs]
  packages_unable_to_install <- setdiff(new_packages, c(needed_cran_pkgs, needed_bioc_pkgs))

  message("Cran packages to install:\n",
           paste0(needed_cran_pkgs, "\n"))

  message("Bioc packages to install:\n",
          paste0(needed_bioc_pkgs, "\n"))

  message("Unknown packages:\n",
          paste0(packages_unable_to_install, "\n"))

  if(do_install){
    if(length(needed_cran_pkgs)) {
      install.packages(needed_cran_pkgs)
    }

    if(length(needed_bioc_pkgs)){
      BiocManager::install(needed_bioc_pkgs)
    }
  }
}
