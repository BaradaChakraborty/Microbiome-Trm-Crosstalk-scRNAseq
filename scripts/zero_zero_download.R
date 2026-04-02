message("Step 1: Pointing R to the local archive file")
dest_file <- "data/raw/GSE205549_RAW.tar"

message("Step 2: Unpacking the archive (R might freeze for a minute, let it work)")
untar(tarfile = dest_file, exdir = "data/raw/")

message("Step 3: Listing the extracted files to confirm success")
list.files("data/raw/")