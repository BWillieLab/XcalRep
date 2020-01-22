#' HR-pQCt phantom measure calibration and reproducibility dataset
#'
#' @source Shriners Hospital for Children, Mereo BioPharma
#' @format A data frame with columns:
#' \describe{
#' \item{site}{site/scanner used to obtain meaures (e.g., Montreal, Toronto, etc.)}
#' \item{timePoint}{time of scan; numeric (e.g., 0, 6, 12)}
#' \item{section}{phantom section; numeric (e.g., value 1-4)}
#' \item{parameter}{measured parameter (e.g., Tb. vBMD, Ct. Ar., etc)}
#' \item{value}{measured value; numeric}
#' \item{phantom}{phantom type (e.g., EFP, QC1)}
#' \item{scanDate}{date of scan (e.g., 2018-02-27)}
#' }
#' @examples
#' \dontrun{
#' # load data into global enviroment
#'  df <- df
#' }
"df"


#' HR-pQCt phantom measure calibration and reproducibility dataset (version 2)
#'
#' @source Shriners Hospital for Children, Mereo BioPharma
#' @format A data frame with columns:
#' \describe{
#' \item{site}{site/scanner used to obtain meaures (e.g., Montreal, Toronto, etc.)}
#' \item{timePoint}{time of scan; numeric (e.g., 0, 6, 12)}
#' \item{scanNumber}{scan ID for given time point for given instrument}
#' \item{section}{phantom section; numeric (e.g., value 1-4)}
#' \item{parameter}{measured parameter (e.g., Tb. vBMD, Ct. Ar., etc)}
#' \item{value}{measured value; numeric}
#' \item{phantom}{phantom type (e.g., EFP, QC1)}
#' \item{scanDate}{date of scan (e.g., 2018-02-27)}
#' }
#' @examples
#' \dontrun{
#' # load data into global enviroment
#'  df <- raw.data
#' }
"raw.data"
