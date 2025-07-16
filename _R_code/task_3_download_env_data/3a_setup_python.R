
#' -----------------------------------------------------------------------------
#' Load `reticulate` package
#' -----------------------------------------------------------------------------
library(reticulate) 

## Run this if Python is not installed on the machine
# install_python() 

reticulate::virtualenv_create(envname = "cm_py")
reticulate::virtualenv_install("cm_py", packages = c("copernicusmarine"))
reticulate::use_virtualenv("cm_py", required = TRUE)
cm <- reticulate::import("copernicusmarine")

## Run once with your username and password #----
cm$login("username","pw")



