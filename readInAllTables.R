library(RMySQL)
library(here)
sisal = dbConnect(MySQL(),user = "root",dbname = 'sisal')#connect to database (only works if you have MySQL setup and the Sisal my sql database)

allFields <- dbListTables(sisal) #get all the fields

all <- vector(mode = "list",length = length(allFields))
for(i in 1:length(all)){#loop through and get all the data
  all[[allFields[i]]] <- dbReadTable(sisal, name = allFields[i])
}

save(list = "all",file = here("allTables.RData"))
