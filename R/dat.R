#' Mayo Clinic Primary Biliary Cirrhosis Data
#' @description The data is from the Mayo Clinic trial in primary biliary cirrhosis (PBC) 
#'    of the liver conducted between 1974 and 1984. A total of 424 PBC patients, referred 
#'    to Mayo Clinic during that ten-year interval, met eligibility criteria for the randomized 
#'    placebo controlled trial of the drug D-penicillamine. The first 312 cases in the data set 
#'    participated in the randomized trial and contain largely complete data. 
#'    This longitudinal data set contains multiple laboratory results, but only on the first 312 patients.
#' @keywords datasets
#' @docType data
#' @usage data(pbc.long)
#' @format A data frame with 1945 rows.
#' \itemize{ 
#' \item \code{id}:  subject identifier
#' \item \code{st}: number of years between registration and the earlier of death, 
#'         transplantion, or study analysis time 
#' \item \code{status}: 0 = alive, 1 = transplanted, 2 = dead 
#' \item \code{status2}: 1 = dead, 0 = alive or transplanted 
#' \item \code{drug}: 1= D-penicillamine, 0=placebo 
#' \item \code{age}: age in years, at registration 
#' \item \code{sex}: 0 = male, 1 = female 
#' \item \code{futime}: number of years between enrollment and this visit date, 
#'         remaining values on the line of data refer to this visit 
#' \item \code{ascites}: presence of ascites: 0 = no, 1 = yes 
#' \item \code{hepatomegaly}: presence of hepatomegaly: 0 = no, 1 = yes 
#' \item \code{spiders}: presence of spiders: 0 = no, 1 = yes 
#' \item \code{edema}: presence of edema: 0 = no edema and no diuretic therapy for edema; 
#'         0.5 = edema present without diuretics, or edema resolved by diuretics; 
#'         1 = edema despite diuretic therapy 
#' \item \code{bilirubin}: serum bilirubin in mg/dl 
#' \item \code{cholesterol}: serum cholesterol in mg/dl 
#' \item \code{albumin}: albumin in gm/dl
#' \item \code{alkaline}: alkaline phosphatase in U/liter
#' \item \code{sgot}: SGOT in U/ml (serum glutamic-oxaloacetic transaminase, 
#'        the enzyme name has subsequently changed to "ALT" in the medical literature) 
#' \item \code{platelets}: platelets per cubic ml/1000 
#' \item \code{proTime}: prothrombin time in seconds 
#' \item \code{histStage}: histologic stage of disease 
#'   }
#' @references Therneau T and Grambsch P, Modeling survival data: extending the Cox Model, Springer-Verlag, 2000.
#' @source \url{https://www.mayo.edu/research/documents/pbcseqhtml/doc-10027141}
#' @examples 
#' summary(pbc.long)
"pbc.long"
########################################################
########################################################
#' Mayo Clinic Primary Biliary Cirrhosis Data
#' @description The data is from the Mayo Clinic trial in primary biliary cirrhosis (PBC) 
#'    of the liver conducted between 1974 and 1984. A total of 424 PBC patients, referred 
#'    to Mayo Clinic during that ten-year interval, met eligibility criteria for the randomized 
#'    placebo controlled trial of the drug D-penicillamine. The first 312 cases in the data set 
#'    participated in the randomized trial and contain largely complete data. 
#'    This data frame contains the first measurement for each patient from the
#'    data set \code{pbc.long}.
#' @keywords datasets
#' @docType data
#' @usage data(pbc.surv)
#' @format A data frame with 312 rows.
#' \itemize{ 
#' \item \code{id}:  subject identifier
#' \item \code{st}: number of years between registration and the earlier of death, 
#'         transplantion, or study analysis time 
#' \item \code{status}: 0 = alive, 1 = transplanted, 2 = dead 
#' \item \code{status2}: 1 = dead, 0 = alive or transplanted 
#' \item \code{drug}: 1= D-penicillamine, 0=placebo 
#' \item \code{age}: age in years, at registration 
#' \item \code{sex}: 0 = male, 1 = female 
#' \item \code{futime}: number of years between enrollment and this visit date, 
#'         remaining values on the line of data refer to this visit 
#' \item \code{ascites}: presence of ascites: 0 = no, 1 = yes 
#' \item \code{hepatomegaly}: presence of hepatomegaly: 0 = no, 1 = yes 
#' \item \code{spiders}: presence of spiders: 0 = no, 1 = yes 
#' \item \code{edema}: presence of edema: 0 = no edema and no diuretic therapy for edema; 
#'         0.5 = edema present without diuretics, or edema resolved by diuretics; 
#'         1 = edema despite diuretic therapy 
#' \item \code{bilirubin}: serum bilirubin in mg/dl 
#' \item \code{cholesterol}: serum cholesterol in mg/dl 
#' \item \code{albumin}: albumin in gm/dl
#' \item \code{alkaline}: alkaline phosphatase in U/liter
#' \item \code{sgot}: SGOT in U/ml (serum glutamic-oxaloacetic transaminase, 
#'        the enzyme name has subsequently changed to "ALT" in the medical literature) 
#' \item \code{platelets}: platelets per cubic ml/1000 
#' \item \code{proTime}: prothrombin time in seconds 
#' \item \code{histStage}: histologic stage of disease 
#'   }
#' @references Therneau T and Grambsch P, Modeling survival data: extending the Cox Model, Springer-Verlag, 2000.
#' @source \url{https://www.mayo.edu/research/documents/pbcseqhtml/doc-10027141}
#' @examples 
#' summary(pbc.surv)
"pbc.surv"