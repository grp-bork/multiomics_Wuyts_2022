# used to install missing packages
if (!require("pacman")) {
    install.packages("pacman", repos = "http://cran.us.r-project.org")
}

# hidden dependency for PDF export, so load it here
p_load("Cairo")
