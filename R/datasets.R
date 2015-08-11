#' Data1
#' 
#' Data1 dataset used to demonstrate the use of SNFtool.
#' 
#' 
#' @name Data1
#' @docType data
#' @format A data frame with 200 observations on the following 2 variables.
#' \describe{ \item{list("V1")}{a numeric vector} \item{list("V2")}{a numeric
#' vector} }
#' @keywords datasets
#' @examples
#' 
#' data(Data1)
#' 
'Data1'

#' Data2
#' 
#' Data2 dataset used to demonstrate the use of SNFtool.
#' 
#' 
#' @name Data2
#' @docType data
#' @format A data frame with 200 observations on the following 2 variables.
#' \describe{ \item{list("V3")}{a numeric vector} \item{list("V4")}{a numeric
#' vector} }
#' @keywords datasets
#' @examples
#' 
#' data(Data2)
#' 
'Data2'

#' dataL
#' 
#' Dataset used to provide an example of predicting the new labels with label
#' propagation.
#' 
#' 
#' @name dataL
#' @docType data
#' @format The format is: List of 2 $ : num [1:600, 1:76] 0.0659 0.0491 0.0342
#' 0.0623 0.062 ...  ..- attr(*, "dimnames")=List of 2 .. ..$ : chr [1:600]
#' "V1" "V2" "V3" "V4" ...  .. ..$ : NULL $ : int [1:600, 1:240] 0 0 0 0 0 0 0
#' 0 0 0 ...  ..- attr(*, "dimnames")=List of 2 .. ..$ : chr [1:600] "V1" "V2"
#' "V3" "V4" ...  .. ..$ : NULL
#' @keywords datasets
#' @examples
#' 
#' data(dataL)
#' 
'dataL'

#' Labels for dataL dataset
#' 
#' The ground truth for dataL dataset
#' 
#' 
#' @name label
#' @docType data
#' @format The format is: int [1:600] 1 1 1 1 1 1 1 1 1 1 ...
#' @keywords datasets
#' @examples
#' 
#' data(label)
#' 
'label'

