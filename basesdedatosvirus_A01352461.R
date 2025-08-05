library(seqinr)
library(stringr)
library(Biostrings)

leer <- function(archivo) {
  secuencia <- readDNAStringSet(archivo)[[1]]
  secuencia_caracter <- as.character(secuencia)
  return(secuencia_caracter)
}

tamaño_secuencia<-function(txt){
  tamaño <- nchar(gsub("[^A-Za-z]", "", txt))
  return(tamaño)
}

cant_nuc <- function(txt, c) {
  
  vec <- gregexpr(c, txt, fixed = TRUE)[[1]]
  cant <- length(vec)
  cantidad <- cant / (tamaño_secuencia(txt)) * 100
  return(cantidad)
}

porc_gc <- function(txt) {
  vec_c <- gregexpr('c', txt, fixed = TRUE)[[1]]
  vec_g <- gregexpr('g', txt, fixed = TRUE)[[1]]
  cant_c <- length(vec_c)
  cant_g <- length(vec_g)
  cant <- cant_g + cant_c
  cantidad <- cant / (tamaño_secuencia(txt)) * 100
  return(cantidad)
}

hebras <- function (txt){
    complementaria <- ""
    for (base in strsplit(txt, NULL)[[1]]) {
      if (base == "A") {
        complementaria <- paste0(complementaria, "T")
      } else if (base == "T") {
        complementaria <- paste0(complementaria, "A")
      } else if (base == "C") {
        complementaria <- paste0(complementaria, "G")
      } else if (base == "G") {
        complementaria <- paste0(complementaria, "C")
      } else {
        complementaria <- paste0(complementaria, base)
      }
    }
    
    return  (complementaria)
  }


dengue <- read.fasta(file = "Dengue.fasta")
mec <- read.fasta(file = "M_corona.fasta")
pneumonia <- read.fasta(file = "Pneu.fasta")
sars <- read.fasta(file = "SARS.fasta")
zika <- read.fasta (file = "zika.fasta")
ha_h1 <- read.fasta (file = "ha_h1n1.fasta")
na_h1 <- read.fasta (file = "na_h1n1.fasta")
ns1_h1 <- read.fasta (file = "ns1_h1n1.fasta")
m1_h1 <- read.fasta (file = "m1_h1n1.fasta")
pa_h1 <- read.fasta (file = "pa_h1n1.fasta")
np_h1 <- read.fasta (file = "np_h1n1.fasta")
pb1_h1 <- read.fasta (file = "pb1_h1n1.fasta")
pb2_h1 <- read.fasta (file = "pb2_h1n1.fasta")

virus <- list (dengue, mec, pneumonia, sars, zika, ha_h1, na_h1, 
               ns1_h1, m1_h1, pa_h1, np_h1, pb1_h1, pb2_h1)
archivos <- list ("Dengue.fasta", "M_corona.fasta", "Pneu.fasta",
                  "SARS.fasta", "zika.fasta", "ha_h1n1.fasta", "na_h1n1.fasta", 
                  "ns1_h1n1.fasta", "m1_h1n1.fasta", "pa_h1n1.fasta", 
                  "np_h1n1.fasta", "pb1_h1n1.fasta", "pb2_h1n1.fasta")

nombres <- c("Dengue", "Middle East respiratory syndrome coronavirus", 
             "Pneumonia", "SARS Coronavirus", "Zika virus", "HA_H1N1", "NA_H1N1",
             "NSI_H1N1", "MI_H1N1", "PA_H1N1", "NP_H1N1", "PB1_H1N1", "PB2_H1N1")

tamaños <- matrix(nrow = 13, ncol = 1)
i <- 1
for (secuencia in virus) {
  tamaños[i,1] <- tamaño_secuencia(secuencia)
  i <- i + 1
}

cant_nuc_a <- matrix(nrow = 13, ncol = 1)
i <- 1
for (secuencia in virus) {
  cant_nuc_a[i, 1] <- cant_nuc(secuencia, 'a')
  i <- i + 1
}

cant_nuc_c <- matrix(nrow = 13, ncol = 1)
i <- 1
for (secuencia in virus) {
  cant_nuc_c[i, 1] <- cant_nuc(secuencia, 'c')
  i <- i + 1
}

cant_nuc_g <- matrix(nrow = 13, ncol = 1)
i <- 1
for (secuencia in virus) {
  cant_nuc_g[i, 1] <- cant_nuc(secuencia, 'g')
  i <- i + 1
}

cant_nuc_t <- matrix(nrow = 13, ncol = 1)
i <- 1
for (secuencia in virus) {
  cant_nuc_t[i, 1] <- cant_nuc(secuencia, 't')
  i <- i + 1
}

p_gc <- matrix(nrow = 13, ncol = 1)
i <- 1
for (secuencia in virus) {
  p_gc[i, 1] <- porc_gc(secuencia)
  i <- i + 1
}

leer_ar <- matrix(nrow = 13, ncol = 1)
compl <- matrix(nrow = 13, ncol = 1)
i <- 1
for (secuencia in archivos){
  leer_ar [i,1] <- substr(leer(secuencia),1,30)
  i <- i + 1
}


i <- 1
for (secuencia in leer_ar) {
  compl[i, 1] <- hebras(secuencia)
  i <- i + 1
}

x <- 1:4
 
dengue_pl <- c(cant_nuc_a [1],cant_nuc_c [1],cant_nuc_g [1],cant_nuc_t [1])
mec_pl <- c(cant_nuc_a [2],cant_nuc_c [2],cant_nuc_g [2],cant_nuc_t [2])
pneumonia_pl <- c(cant_nuc_a [3],cant_nuc_c [3],cant_nuc_g [3],cant_nuc_t [3])
sars_pl <- c(cant_nuc_a [4],cant_nuc_c [4],cant_nuc_g [4],cant_nuc_t [4])
zika_pl <- c(cant_nuc_a [5],cant_nuc_c [5],cant_nuc_g [5],cant_nuc_t [5])
ha_h1_pl <- c(cant_nuc_a [6],cant_nuc_c [6],cant_nuc_g [6],cant_nuc_t [6])
na_h1_pl <- c(cant_nuc_a [7],cant_nuc_c [7],cant_nuc_g [7],cant_nuc_t [7])
ns1_h1_pl <- c(cant_nuc_a [8],cant_nuc_c [8],cant_nuc_g [8],cant_nuc_t [8])
m1_h1_pl <- c(cant_nuc_a [9],cant_nuc_c [9],cant_nuc_g [9],cant_nuc_t [9])
pa_h1_pl <- c(cant_nuc_a [10],cant_nuc_c [10],cant_nuc_g [10],cant_nuc_t [10])
np_h1_pl <- c(cant_nuc_a [11],cant_nuc_c [11],cant_nuc_g [11],cant_nuc_t [11])
pb1_h1_pl <- c(cant_nuc_a [12],cant_nuc_c [12],cant_nuc_g [12],cant_nuc_t [12])
pb2_h1_pl <- c(cant_nuc_a [13],cant_nuc_c [13],cant_nuc_g [13],cant_nuc_t [13])


plot ("Nucleótido","Cantidad de nulceótido", xlim = c(1,4) , ylim = c(0,40) , 
      type = "l", main = "composiciones de nucleótidos")

lines (x, dengue_pl, col = "blue")
lines (x, mec_pl, col = "purple")
lines (x, pneumonia_pl, col = "red")
lines (x, sars_pl, col = "yellow")
lines (x, zika_pl, col = "green")
lines (x, ha_h1_pl, col = "brown")
lines (x, na_h1_pl, col = "gray")

legend("bottomleft", legend = c("Dengue", "Middle East respiratory syndrome coronavirus",
                              "Pneumonia", "SARS Coronavirus", "Zika virus", "HA_H1N1", 
                              "NA_H1N1"), col = c("blue", "purple", "red", "yellow", 
                                                  "green", "brown", "gray"), lwd = 2, bty = "n")

plot ("Nucleótido","Cantidad de nulceótido", xlim = c(1,4) , ylim = c(0,40) , 
      type = "l", main = "composiciones de nucleótidos")

lines (x, ns1_h1_pl, col = "blue")
lines (x, m1_h1_pl, col = "purple")
lines (x, pa_h1_pl, col = "red")
lines (x, np_h1_pl, col = "yellow")
lines (x, pb1_h1_pl, col = "green")
lines (x, pb2_h1_pl, col = "brown")

legend("bottomleft", legend = c("NS1_H1N1", "M1_H1N1","PA_H1N1", "NP_H1N1", "PB1_H1N1", "PB2_H1N1", 
                              "NA_H1N1"), col = c("blue", "purple", "red", "yellow", 
                                                  "green", "brown"), lwd = 2, bty = "n")

datos_tama <- data.frame("Nombre" = nombres,"Tamaño de las secuencias." = tamaños, "Porcentaje de la adenina" = cant_nuc_a, "Porcentaje de la timina" = cant_nuc_t, "Porcentaje de la guanina" = cant_nuc_g, "Porcentaje de la citocina" = cant_nuc_c, "Porcentaje de GC" = p_gc, "Hebra directa" = leer_ar, "Hebra complementaria" = compl)
print (datos_tama)
