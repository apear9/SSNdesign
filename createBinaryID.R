createBinaryID <- function (ssn, o.write) 
{
  if (file.exists("binaryID.db") == TRUE) {
    if (o.write == TRUE) {
      unlink("binaryID.db")
    }
    else {
      cat("binaryID.db already exists - no changes were made to binaryID.db table\n")
      mm <- T
    }
  }
  options(show.error.messages = FALSE)
  m <- try(if (file.exists("binaryID.db") == FALSE) {
    mm <- F
    driver <- RSQLite::SQLite()
    db.name <- "binaryID.db"
    connect <- dbConnect(SQLite(), db.name)
    net.no <- as.numeric(levels(ssn@network.line.coords$NetworkID))
    for (i in 1:length(net.no)) {
      network <- paste("net", net.no[i], sep = "")
      file.name <- paste("netID", net.no[i], ".dat", sep = "")
      if (dbExistsTable(connect, network)) {
        dbRemoveTable(connect, network)
      }
      dbWriteTable(connect, network, read.table(file = file.name, 
                                                header = T, sep = ",", colClasses = c("numeric", 
                                                                                      "character")), overwrite = T, row.names = F)
      if (i == length(net.no)) {
        if (length(dbListTables(connect)) != length(net.no)) {
          dbDisconnect(connect)
          stop("ERROR: binary tables did not import to SQLite database properly")
        }
      }
    }
    dbDisconnect(connect)
  }, silent = TRUE)
  options(show.error.messages = TRUE)
  if (mm != T) {
    if (m != T) {
      dbDisconnect(connect)
      stop("ERROR: binary tables did not import to SQLite database properly")
    }
  }
}