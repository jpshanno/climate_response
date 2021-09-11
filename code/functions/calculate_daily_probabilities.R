##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param in.path
##' @return
##' @author Joe Shannon
##' @export
calculate_daily_probabilities <- function(in.path) {

  cmd <- 
    glue::glue(
      "zcat <in.path> | 
      awk \\
          -F',' ' 
          BEGIN{
            SUBSEP=\",\" 
            print(\"site,site_status,gcm,scenario,simulation_date,prop_above_0,prop_above_neg_10,prop_above_neg_25,prop_above_neg_50,prop_above_neg_100\")
          } 
          NR>1{
            prop_0[$2,$3,$4,$5,$6]+=($7>0)
            prop_10[$2,$3,$4,$5,$6]+=($7>-10)
            prop_25[$2,$3,$4,$5,$6]+=($7>-25)
            prop_50[$2,$3,$4,$5,$6]+=($7>-50)
            prop_100[$2,$3,$4,$5,$6]+=($7>-100)
            n[$2,$3,$4,$5,$6]++
          } 
          END{ 
            for(i in prop_10) print i\",\"prop_0[i]/n[i]\",\"prop_10[i]/n[i]\",\"prop_25[i]/n[i]\",\"prop_50[i]/n[i]\",\"prop_100[i]/n[i] | \"sort\"
          } '",
      .open = "<",
      .close = ">")
  
  fread(cmd = cmd,
        colClasses = list(character = "site", 
                          Date = "simulation_date"),
        nThread = 4)

}
