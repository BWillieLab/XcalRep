#' HR-pQCT QC1 phantom calibration and reproducibility dataset
#'
#' @source Shriners Hospital for Children, Mereo BioPharma
#' @format A data frame with columns:
#' \describe{
#' \item{scanID}{Scan identfier}
#' \item{site}{Site/scanner used to obtain meaures (e.g., Montreal, Toronto, etc.)}
#' \item{parameter}{Scanned parameter (e.g., Tb. vBMD, Ct. Ar., etc)}
#' \item{timePoint}{Time of scan; numeric (e.g., 0, 6, 12)}
#' \item{value}{Measured value; numeric}
#' \item{section}{Phantom section; numeric (e.g., value 1-4)}
#' \item{scanner}{HR-pQCT scanner (e.g., XCT, XCT2)}
#' }
#' @examples
#' \dontrun{
#' # load data into global enviroment
#'  df <- qc1
#' }
"qc1"


#' HR-pQCT EFP phantom calibration and reproducibility dataset
#'
#' @source Shriners Hospital for Children, Mereo BioPharma
#' @format A data frame with columns:
#' \describe{
#' \item{scanID}{Scan identfier}
#' \item{site}{Site/scanner used to obtain meaures (e.g., Montreal, Toronto, etc.)}
#' \item{parameter}{Scanned parameter (e.g., Tb. vBMD, Ct. Ar., etc)}
#' \item{timePoint}{Time of scan; numeric (e.g., 0, 6, 12)}
#' \item{value}{Measured value; numeric}
#' \item{section}{Phantom section; numeric (e.g., value 1-4)}
#' \item{scanner}{HR-pQCT scanner (e.g., XCT, XCT2)}
#' }
#' @examples
#' \dontrun{
#' # load data into global enviroment
#'  df <- efp
#' }
"efp"

#' HR-pQCT in vivo patient reproducibility dataset
#'
#' @source Shriners Hospital for Children, Mereo BioPharma
#' @format A data frame with columns:
#' \describe{
#' \item{scanID}{Scan identfier}
#' \item{site}{Site/scanner used to obtain meaures (e.g., Montreal, Toronto, etc.)}
#' \item{parameter}{Scanned parameter (e.g., Tb. vBMD, Ct. Ar., etc)}
#' \item{timePoint}{Time of scan; numeric (e.g., 0, 6, 12)}
#' \item{value}{Measured value; numeric}
#' \item{scanner}{HR-pQCT scanner (e.g., XCT, XCT2)}
#' }
#' @examples
#' \dontrun{
#' # load data into global enviroment
#'  df <- invivo
#' }
"invivo"

#' HR-pQCT in vivo patient reproducibility dataset (unregistered)
#'
#' @source Shriners Hospital for Children, Mereo BioPharma
#' @format A data frame with columns:
#' \describe{
#' \item{scanID}{Scan identfier}
#' \item{site}{Site/scanner used to obtain meaures (e.g., Montreal, Toronto, etc.)}
#' \item{parameter}{Scanned parameter (e.g., Tb. vBMD, Ct. Ar., etc)}
#' \item{timePoint}{Time of scan; numeric (e.g., 0, 6, 12)}
#' \item{value}{Measured value; numeric}
#' \item{scanner}{HR-pQCT scanner (e.g., XCT, XCT2)}
#' \item{replicateSet}
#' \item{age}{Patient age (years); numeric}
#' \item{sex}{Patient sex; character (e.g., Male, Female)}
#' \item{weight}{Patient weight (kg); numeric}
#' \item{registration}{Registration status; character (e.g., CAS.registered, unregistered)}
#' }
#' @examples
#' \dontrun{
#' # load data into global enviroment
#'  df <- invivo_unregistered
#' }
"invivo_unregistered"

#' HR-pQCT in vivo patient reproducibility dataset (registered)
#'
#' @source Shriners Hospital for Children, Mereo BioPharma
#' @format A data frame with columns:
#' \describe{
#' \item{scanID}{Scan identfier}
#' \item{site}{Site/scanner used to obtain meaures (e.g., Montreal, Toronto, etc.)}
#' \item{parameter}{Scanned parameter (e.g., Tb. vBMD, Ct. Ar., etc)}
#' \item{timePoint}{Time of scan; numeric (e.g., 0, 6, 12)}
#' \item{value}{Measured value; numeric}
#' \item{scanner}{HR-pQCT scanner (e.g., XCT, XCT2)}
#' \item{replicateSet}
#' \item{age}{Patient age (years); numeric}
#' \item{sex}{Patient sex; character (e.g., Male, Female)}
#' \item{weight}{Patient weight (kg); numeric}
#' \item{registration}{Registration status; character (e.g., CAS.registered, unregistered)}
#' }
#' @examples
#' \dontrun{
#' # load data into global enviroment
#'  df <- invivo_registered
#' }
"invivo_registered"
