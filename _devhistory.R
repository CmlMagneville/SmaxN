usethis::use_build_ignore("_devhistory.R")

usethis::use_git(message = ":tada: Initial commit")

usethis::edit_file("DESCRIPTION")

usethis::use_git(message = ":bulb: Edit package metadata")

usethis::use_package_doc()