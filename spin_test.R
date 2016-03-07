#+ results='asis', echo = FALSE
for (i in 1:5) {
  cat("## This is a heading for ", i, "\n")
  cat("<!-- This is a comment for ", i, "-->\n")
  print(i)    
}