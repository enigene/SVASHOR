#!/usr/bin/env Rscript
# Part of svashor which generate reports with Structural Variants of
# Alpha Satellite Higher-Order Repeats.
# Inputs *-HOR.tsv with structural variants, *-stat.txt with counted
# variants, and table ASMonColors.txt with color information.
argv <- commandArgs(T)

library("xtable")
library("xml2")
library("stringr")

vtab <- read.delim(argv[1], header = T, stringsAsFactors = T, strip.white = T)
stab <- read.delim(argv[2], header = F, stringsAsFactors = T, strip.white = T)
outf <- file(argv[3], "wa")

# write css style
writeLines("<link href=\"style.css\" rel=\"stylesheet\" media=\"all\">", outf)

# Replace M1+ to M1 in column Class
levels(vtab$Class)[levels(vtab$Class)=="M1+"] <- "M1"

# Variant tables
facHVar <- with(vtab, factor(HVar, levels = unique(HVar)))
facPVar <- with(vtab, factor(PVar, levels = unique(PVar)))

g1 <- split(vtab, with(vtab, interaction(facHVar,facPVar)), drop = T)
g2 <- lapply(g1, t)

for(i in 1:length(g2)){
  frag = g2[[i]][1]
  HORAbbr = g2[[i]][2]
  Hvar = g2[[i]][3]
  Pvar = g2[[i]][4]
  hv = sub("var", "", Hvar)
  pv = sub("var", "", Pvar)
  vtabCaption = paste0("Variant #", hv, " of HOR ", HORAbbr, " in ", frag)
  sub.g2 = as.data.frame(g2[[i]][-c(1:4),])
  names(sub.g2) <- NULL

  subtable <- print.xtable(xtable(sub.g2),
              html.table.attributes = NULL,
              type = "html",
              comment = F,
              print.results = F,
              include.colnames = F)

 subtable.xml <- read_xml(subtable)
 xml_add_sibling(
   xml_find_first(subtable.xml, "//table"),
   .value = as_xml_document(paste0("<span class=\"tblDesc\">", vtabCaption, "</span>")),
   .where = "before"
 )
 xml_add_sibling(
   xml_find_first(subtable.xml, "//table"),
   .value = as_xml_document(paste0("<br />")),
   .where = "after"
 )
 allNodes <- xml_find_all(subtable.xml, "//td")
 lapply(allNodes, function(x) {
   xml_set_text(x, xml_text(x, trim = T))
 })
 classNodes <- xml_find_all(subtable.xml, "//td[text()='Class']/following-sibling::td")
 lapply(classNodes, function(x) {
   xml_set_attr(x, "class", tolower(xml_text(x, trim = T)))
 })
 classNodes <- xml_find_all(subtable.xml, "//table")
 lapply(classNodes, function(x) {
   xml_set_attr(x, "class", "tblVar")
 })
 write_xml(subtable.xml, outf, options = c("no_declaration"))
}


# Table with structural variant statistics

# Omit rows if in the first column is present single value
stab <- stab[(str_split_fixed(stab$V1," " ,2)[,2])!="", ]
rownames(stab) <- NULL

# Stat table constants:
# threshold for variants
thld <- 5
# limiting number of rows
rowsLim <- 10

# vector with counted values
HORcnt <- stab$V2

topItem <- HORcnt[1]
numItems <- length(HORcnt)
thldRow <- length(which(HORcnt>=thld))

if((thldRow+1) < rowsLim) {
  stabs <- na.omit(rbind(stab[HORcnt>=thld, ], stab[c((thldRow+1):rowsLim), ]))
} else {
  stabs <- stab[HORcnt>=thld, ]
}

stabs <- cbind(Var=rownames(stabs), stabs)

stabsCaption = paste0("Structural variants of HOR ", HORAbbr, " in ", frag)

colnames(stabs) <- c("Var#","Sequence of monomers in a HOR","Number of copies")

stabsx <- print.xtable(
  xtable(stabs),
  html.table.attributes = NULL,
  type = "html",
  comment = F,
  print.results = F,
  include.rownames = F
)

stabsx.xml <- read_xml(stabsx)

xml_add_sibling(
  xml_find_first(stabsx.xml, "//table"),
  .value = as_xml_document(paste0("<span class=\"tblDesc\">", stabsCaption, "</span>")),
  .where = "before"
)
allNodes <- xml_find_all(stabsx.xml, "//td")
invisible(lapply(allNodes, function(x) {
  xml_set_text(x, xml_text(x, trim = T))
}))
classNodes <- xml_find_all(stabsx.xml, "//table")
invisible(lapply(classNodes, function(x) {
  xml_set_attr(x, "class", "tblStatr")
}))
xml_add_sibling(
  xml_find_first(stabsx.xml, "//table"),
  .value = as_xml_document(paste0("<hr class=\"tblHr\" />")),
  .where = "after"
)

write_xml(stabsx.xml, outf, options = c("no_declaration"))

close(outf)
