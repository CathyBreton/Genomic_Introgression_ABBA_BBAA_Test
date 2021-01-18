
#################################################################################################################
###   Format_ABBA_Table_Musa  ###
###   Author : Catherine Breton ###
###   Name :
###   Date : 01/29/2020 11:36 am
###   Version : 
#     Created : 12/20/2018 08:00 am
#     Revision : 
###   R version developpement : 3.62
### 
###
#     Intellectual property belongs to Bioversity
#     Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html)
#################################################################################################################

# Remove R memory
rm(list=ls())

#################################################################################################################
##### Data Preparation
#################################################################################################################

# Load and install package 


####################################################


##### Data Formated Table 



###################################################





# Get the command-line arguments.
args <- commandArgs(trailingOnly = TRUE)
allele_freqs_file_name <- args[1]
ABBA_table_input  <- args[1]
output_file_name <- args[2]
spc_p1 <- args[3]
spc_p2 <- args[4]
spc_p3 <- args[5]
spc_o <- args[6]
lg_length_table_file_name <- args[7]
output2_file_name <- args[8]




### Instal Formattable package
install.packages("formattable") 
# Charger la librarie avant de l'utiliser 
library(formattable) 
# Charger les données

## Load Data
setwd("~/Bioversity_Work/2020_ABBA_Julie/ABBA_Python_Martin_Banana/ABBA_Python_Result")

# Read the Data Table 

Data <- read.table(file.choose(), header = FALSE, sep = "")
View(Data)
ABBA_table <- Data
#ABBA_table = read.table(ABBA_table_input)
ABBA_table = read.table("P3A_TAB_abba_baba.txt")
# Remove line empty
ABBA_table[!apply(ABBA_table == "", 1, all),]
# Remove the first colonne 
ABBA_table <- ABBA_table[,-1] # On supprime la 1eme colonne

names(ABBA_table) <- c("Pop1", "Pop2", "pop3", "Outgroup", "D", "D_sd", "D_err", "D_Z", "D_p", "fd", "ABBA", "BABA", "BBAA")

ABBA_table <- data.frame(ABBA_table)
formattable(df, list(area(col = a:c) ~ color_tile("transparent", "pink")))


# Sélectionner quelques colonnes comme dans le tableau Excel 
tab_reduce = ABBA_table[, c(1, 2, 3, 4, 5, 8)] 

# Créer un tableau formattable 
widget_formattable = formattable(tab_reduce) 
# Afficher dans un navigateur web ou dans RStudio (si vous l'utilisez) 
widget_formattable 

###################################################################################

##Correction example 

# Un message d'erreur nous avertit qu'il y a des valeurs 'None' et qu'elles seront remplacées par des NA. -> C'est correct.
tab_reduce$gene_FPKM = as.numeric(tab_reduce$gene_FPKM) 
# Utiliser 'accounting' pour formater les valeurs numériques  
tab_reduce$gene_FPKM = accounting(tab_reduce$gene_FPKM) 




# Utilisation de fonctions 'built-in' pour afficher des élements graphiques. 
widget_formattable = formattable(tab_reduce, list(
  D = color_tile('lightblue', 'white'),
  # Nous avons besoin de mettre na.rm=TRUE pour gérer les valeurs NA présentes dans les données
  D_Z = color_bar('red', na.rm = TRUE)
  #mascot_score = color_tile('white', 'orange'),
  #area(col=c(intensity_1, intensity_2, intensity_3, intensity_4)) ~ normalize_bar("pink")
))
widget_formattable


#####################################################################################

install.packages('htmltools')
install.packages('webshot')
library(htmltools)
library(webshot)
webshot::install_phantomjs()

# Exporter le tableau en PNG, PDF ou JPEG
# Voir: https://github.com/renkun-ken/formattable/issues/26
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}





export_formattable(widget_formattable, "./CombinaisainABCDEF_Tableau.png")
# Pour le PDF, les couleurs ne sont pas disponibles.  C'est un problème connu. 
# Voir: https://github.com/renkun-ken/formattable/issues/53
export_formattable(widget_formattable, "./CombinaisainABCDEF_Tableau.pdf")
export_formattable(widget_formattable, "./CombinaisainABCDEF_Tableau.jpeg")



########################################################################################


#####################################
#### Function format the BBA Table
#####################################
formatABBATab <- function (ABBA_Data)
{
  ABBA_table <- ABBA_Data
  # Remove line empty
  ABBA_table[!apply(ABBA_table == "", 1, all),]
  # Remove the first colonne 
  Result <- ABBA_table[,-1] # On supprime la 1eme colonne
  # Remove the last colonne 
  Result <- Result[,-14] # On supprime la 1eme colonne
  #view(Result)
  # Raname the vector
  names(Result) <- c("Pop1", "Pop2", "pop3", "Outgroup", "D", "D_sd", "D_err", "D_Z", "D_p", "fd", "ABBA", "BABA", "BBAA")
  # Transform in dataframe
  Result <- data.frame(Result)
  ## Table
  return(Result)
}

Result <- formatABBATab(ABBA_Data)
View(Result)

write.table(x = Result, file = "366_Combinaison_TAB_abba_baba_PNG_Corrected.csv")

write.table(x = Result, file = "ITC0568_TAB_abba_baba.csv")




