unlink("docs", recursive = TRUE)
bookdown::render_book("docs")
system("touch docs/.nojekyll")
