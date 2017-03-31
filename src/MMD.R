## Source: http://www.anthropology.uw.edu.pl/05/MMD.R

## Author: Arkadiusz Soltysiak

## This script is a supplementary file to the following paper:
## Soltysiak A. (2011), Technical Note: An R script for Smith's Mean Measure of Divergence,
## "Bioarchaeology of the Near East" 5:41-44.

## Paper can be found in the /docs/ folder for this repository.

## -------------SECTION A: INPUT-----------------------------------------------

file_n = choose.files(default = "", caption = "Select a CSV file with numbers of observations", multi = FALSE)
if (file_n == "") stop("The file with numbers of observation was not loaded!")
M_n = read.csv(file=file_n)
if (file_n == "") file_n = "ERROR"

file_p = choose.files(default = "", caption = "Select a CSV file with trait proportions", multi = FALSE)
if (file_p == "") stop("The file with trait proportions was not loaded!")
M_p = read.csv(file=file_p)
if (file_p == "") file_p = "ERROR"

## -------------SECTION B: MMD FORMULA-----------------------------------------

## Anscombe transformation
theta <- function(n,p) { asin((n/(n+3/4))*(1-2*p)) }
## Freeman & Tukey transformation
#theta <- function(n,p) { 0.5*(asin(1-(2*p*n/(n+1)))+asin(1-(2*((p*n)+1)/(n+1)))) }

## Freeman & Tukey correction
thetadiff <- function(nA,pA,nB,pB) { (theta(nA,pA) - theta(nB,pB))^2 - (1/(nA+0.5) + 1/(nB+0.5)) }
## Grewal correction
#thetadiff <- function(nA,pA,nB,pB) { (theta(nA,pA) - theta(nB,pB))^2 - (1/nA + 1/nB) }
## uncorrected formula
#thetadiff <- function(nA,pA,nB,pB) { (theta(nA,pA) - theta(nB,pB))^2 }

## -------------SECTION C: VARIABLE STATUS-------------------------------------

VarMatrix <- M_n[1:2,2:length(M_n[1,])]
MMDMatrix <- matrix(0,length(M_n[,1]),length(M_n[,1]))

for (a in seq_along(VarMatrix[1,])) {
for (b in seq_along(MMDMatrix[,1])) {
for (c in seq_along(MMDMatrix[1,])) {
## percent frequencies
MMDMatrix[b,c] <- thetadiff(M_n[b,a+1], M_p[b,a+1]/100, M_n[c,a+1], M_p[c,a+1]/100)
## proportions
#MMDMatrix[b,c] <- thetadiff(M_n[b,a+1], M_p[b,a+1], M_n[c,a+1], M_p[c,a+1])
} }

for (b in seq_along(MMDMatrix[,1])) {
for (c in seq_along(MMDMatrix[1,])) {
if (b >= c) { MMDMatrix[b,c] = 0 } } }

VNeg <- 0
VPos <- 0
for (b in seq_along(MMDMatrix[,1])) {
for (c in seq_along(MMDMatrix[1,])) {
if (MMDMatrix[b,c] > 0) { VPos = VPos + 1 }
if (MMDMatrix[b,c] < 0) { VNeg = VNeg + 1 } } }

VarMatrix[1,a] = sum(MMDMatrix)
VarMatrix[2,a] = VPos / (VPos + VNeg) }

VarStatus <- t(VarMatrix)

## -------------SECTION D: MMD MATRIX------------------------------------------

MMDMatrix <- matrix(0,length(M_n[,1]),length(M_n[,1]))
dimnames(MMDMatrix) <- list(M_n[,1], M_n[,1])

for (a in seq_along(MMDMatrix[,1])) {
for (b in seq_along(MMDMatrix[1,])) {
MMDVect <- vector("double",length(M_n[1,])-1)
## percent frequencies
for (i in seq_along(MMDVect)) { MMDVect[i] <- thetadiff(M_n[a,i+1], M_p[a,i+1]/100, M_n[b,i+1], M_p[b,i+1]/100) }
## proportions
#for (i in seq_along(MMDVect)) { MMDVect[i] <- thetadiff(M_n[a,i+1], M_p[a,i+1], M_n[b,i+1], M_p[b,i+1]) }
MMDMatrix[a,b] <- sum(MMDVect) / length(MMDVect) } }

## forced 0 when a sample is compared to itself
#for (a in seq_along(MMDMatrix[,1])) { MMDMatrix[a,a] = 0 }

## -------------SECTION E: SD MATRIX-------------------------------------------

## standard deviation for MMD, Sjovold's formula
SDMatrix <- matrix(0,length(M_n[,1]),length(M_n[,1]))
dimnames(SDMatrix) <- list(M_n[,1], M_n[,1])
SDDiff <- function(nA,nB) { (1/nA + 1/nB)^2 }

for (a in seq_along(MMDMatrix[,1])) {
for (b in seq_along(MMDMatrix[1,])) {
SDVect <- vector("double",length(M_n[1,])-1)
for (i in seq_along(SDVect)) { SDVect[i] <- SDDiff(M_n[a,i+1], M_n[b,i+1]) }
SDMatrix[a,b] <- sqrt(sum(SDVect)*2 / length(SDVect)^2) } }

## -------------SECTION F: SIGNIFICANCE MATRIX---------------------------------

## statistical significance
SigMatrix <- matrix(1,length(M_n[,1]),length(M_n[,1]))
dimnames(SigMatrix) <- list(M_n[,1], M_n[,1])

for (a in seq_along(MMDMatrix[,1])) {
for (b in seq_along(MMDMatrix[1,])) {

dist <- MMDMatrix[a,b] / SDMatrix[a,b]
SigMatrix[a,b] = round((1-pnorm(dist))*2, digits = 8)
if (MMDMatrix[a,b] < 0) SigMatrix[a,b] = 1
} }


## -------------SECTION G: OUTPUT----------------------------------------------

VarStatus
#write.csv(VarStatus, file="output_vars.csv")

MMDMatrix
#write.csv(MMDMatrix, file="output_mmd.csv")

SDMatrix
#write.csv(SDMatrix, file="output_mmd-sd.csv")

SigMatrix
#write.csv(SigMatrix, file="output_mmd-p.csv")

cat("The script was executed for files", file_n, "and", file_p, ".", fill = TRUE)