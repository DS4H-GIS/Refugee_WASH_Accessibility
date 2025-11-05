#!/usr/local/bin/Rscript

renv::restore(lockfile="renv.lock", prompt=FALSE)

library(tidyverse)
#library(data.table)
library(sf)
library(terra)
library(tidyterra)
source("./utils/spatial.R")
source("./utils/ACC.R")


get_ACC_D <- function(acc25, acc22,aoigrid){
    
    selected_vars <- c("LT_t", "P_t", "S_t", "AVG_t")
    pop_vars22 <- c("Total22Feb", "Male22Feb", "Female22Feb")
    pop_vars25 <- c("Total25Jan", "Male25Jan", "Female25Jan")
    
    name22 <- c("g_index", pop_vars22, paste0(selected_vars,"_22"))
    name25 <- c("g_index", pop_vars25, paste0(selected_vars,"_25"))
    
    ACC22_tib <- as_tibble(acc22) |>
        select(g_index, all_of(pop_vars22), all_of(selected_vars)) |> 
        rename_at(vars(c("g_index", pop_vars22, selected_vars)),~name22)
    ACC25_tib <- as_tibble(acc25) |> 
        select(g_index, all_of(pop_vars25), all_of(selected_vars)) |> 
        rename_at(vars(c("g_index", pop_vars25, selected_vars)),~name25)
    
    tidyterra::left_join(aoigrid,ACC22_tib) |> tidyterra::left_join(x = _,ACC25_tib,by = join_by(g_index)) |>
        drop_na()|>
        #mutate(across(everything(), ~replace_na(.x, 0))) |>
        mutate(d_pop_t = Total25Jan - Total22Feb,
               d_pop_m = Male25Jan - Male22Feb,
               d_pop_f = Female25Jan - Female22Feb,
               d_LT_t = LT_t_25 - LT_t_22,
               d_P_t = P_t_25 - P_t_22,
               d_S_t = S_t_25 - S_t_22,
               d_AVG_t = AVG_t_25 - AVG_t_22) |> 
        select(d_pop_t, d_pop_m, d_pop_f, d_LT_t, d_P_t, d_S_t, d_AVG_t)
    
}


get_2024_latr <- function(){
    
    library(pdftools)
    
    x <- pdf_text("https://rohingyaresponse.org/wp-content/uploads/2024/11/Overview-and-Monitoring-of-WASH-Per-Camp_Round_5_October-31_2024.pdf")|>
        read_lines()
    
    camp_lines <- x[str_detect(x, "^\\s*Camp\\s")]
    
    parse_camp <- function(s){
        s0 <- str_trim(s)                              # 선·후공백만 제거. 내부 공백 유지
        SMSDCamp <- str_match(s0,
                              "^(Camp\\s+(?:\\d{2}(?:E|W)?(?:\\s+Ext)?|KRC|NRC))")[,1]
        rest <- str_remove(s0, fixed(SMSDCamp))
        vals <- str_split(rest, "\\s{2,}", simplify = TRUE) |> as.vector()
        #stopifnot(length(vals) >= 17)                  # 기대 열수 확인
        
        tibble(
            SMSDCamp = SMSDCamp,
            fstp_cod_ok_pct              = vals[2],
            people_per_functional_latrine= vals[3],
            facility_functionality_pct   = vals[4],
            feel_safe_pct                = vals[5],
            avg_5star_rate               = vals[6],
            e_coli_zero_pct              = vals[7],
            liters_collected_pp_pd       = vals[8],
            tapstand_primary_pct         = vals[9],
            enough_water_pct             = vals[10],
            hhw_collection_gppd          = vals[11],
            visible_waste_pct            = vals[12],
            waste_segregation_pct        = vals[13],
            soap_present_pct             = vals[14],
            mhm_access_pct               = vals[15],
            awd_cases_addressed_pct      = vals[16],
            wash_corrective_actions_pct  = vals[17],
            priority_index               = vals[18]
        )
    }
    
    ret <- map_dfr(camp_lines, parse_camp) |>
        mutate(
            SMSDCamp = str_replace(SMSDCamp, "(?i)^\\s*camp\\s+(\\d{2})\\s+ext\\b", "Camp \\1X"),
            SMSDCamp = str_replace(SMSDCamp, "(?i)^\\s*camp\\b", "Camp")
        )
    
    
    return(ret)
    
}




main <- function(){
    
    #constants
    cellsize <- 50L
    
    cat("\n\n\nThe output files of running this script are located at './out/*' \n\n\n")
    aoi_crs <- rast("./data/img_sample/115258_198184.png") |> terra::crs()
    
    
    cat("\n\n\nReading data from files\n\n\n")
    
    # Make grid out of camp's boundary
    bndry <- vectNreproject("data/result/Camp_100m_buffer.shp",aoi_crs)
    
    
    A1pop <- tryCatch({
        
        vectNreproject("./data/camp_outline/20230412_a1_camp_outlines.kml",aoi_crs) |>
            left_join(read_csv("./data/camp_outline/Population.csv")) |>
            select("CampLabel","SMSDCamp",contains(c("22","25")))
        
    },error = function(e){
        # due to missing LIBKML driver on Windows's GDAL
        vectNreproject("./data/camp_outline/20230412_a1_camp_outlines.gpkg",aoi_crs) |>
            left_join(read_csv("./data/camp_outline/Population.csv")) |>
            select("CampLabel","SMSDCamp",contains(c("22","25")))
        
    })
    
    R2022 <- vectNreproject("/vsizip/./data/result/Rohingya_z18_45441_year2022_2025v7.zip/Rohingya_z18_45441_year2022_v7.gpkg",aoi_crs) |> mutate(i_index = row_number())
    R2025 <- vectNreproject("/vsizip/./data/result/Rohingya_z18_45441_year2022_2025v7.zip/Rohingya_z18_00000_year2025_v7.gpkg",aoi_crs) |> mutate(i_index = row_number())

    
    roads <- st_read("./data/road/20250910_Access_Road_Footpath_all_camps.shp") |> st_transform(aoi_crs)
    
    
    
    
    
    # latrine vsizip/
    f_latr22 <- vectNreproject("/vsizip/./data/facility/Rohingya_refugee_response.zip/Rohingya_refugee_response/WASH_Latrine_20220531.shp",
                               aoi_crs) |>
        _[bndry]|>
        mutate(LT_Male = ifelse(is.nan(LT_Male ), 0, LT_Male)) |>
        mutate(LT_all_gen = ifelse(is.nan(LT_all_gen ), 0, LT_all_gen)) |>
        mutate(LT_Male_sum = LT_all_gen+LT_Male) |>
        mutate(LT_Female = ifelse(is.nan(LT_Female), 0, LT_Female)) |>
        mutate(LT_Female_sum = LT_all_gen+LT_Female)
    
    f_latr24 <- vectNreproject("/vsizip/./data/facility/Rohingya_refugee_response.zip/Rohingya_refugee_response/WASH_Latrine_20240815.shp",aoi_crs) |>
        _[bndry] |>
        mutate(nb_Latrine = as.numeric(nb_Latrine)) 
    
    # Showers
    f_shower22 <- vectNreproject("/vsizip/./data/facility/Rohingya_refugee_response.zip/Rohingya_refugee_response/WASH_Bath_20220531.shp",aoi_crs)|>
        _[bndry] |>
        mutate(Bathing_M = Bathing_M + Bath_gen_u) |> 
        mutate(Bathing_F = Bathing_F + Bath_gen_u)
    f_shower24 <- vectNreproject("/vsizip/./data/facility/Rohingya_refugee_response.zip/Rohingya_refugee_response/WASH_Bath_20240815.shp",aoi_crs)|>
        _[bndry]
    
    # water pumps
    f_pumps22 <- vectNreproject("/vsizip/./data/facility/Rohingya_refugee_response.zip/Rohingya_refugee_response/WASH_handpump_20220531.shp",aoi_crs) |>
        _[bndry] |> mutate(TW = as.numeric(TW))
    f_pumps24 <- vectNreproject("/vsizip/./data/facility/Rohingya_refugee_response.zip/Rohingya_refugee_response/WASH_handpump_20240815.shp",aoi_crs) |>
        _[bndry] |>
        mutate(nb_TW = as.numeric(nb_TW))
    
    cat("Making Grid & Performing Areal interpolation\n\n\n")
    aoigrid <- vect_make_grid(bndry,cellsize = c(cellsize,cellsize),clip = T)
    
    aux_attr <- c("Total22Feb", "Female22Feb", "Male22Feb", "Total25Jan", "Female25Jan", "Male25Jan")
    R2022_grid <- weighted_sum(R2022,aoigrid,aux_vect = A1pop, aux_attr = aux_attr[1:3])
    R2025_grid <- weighted_sum(R2025,aoigrid,aux_vect = A1pop, aux_attr = aux_attr[4:6])
    
    
    
    ############################
    
    # Scenario I
    cat("Calculating Scenario I ACC (1/3)\n\n\n")
    
    cat("Calculating road network distance. Takes some time (Approx. 10 min).\nPlease be patient.\n\n\n")
    Dmat_latr24 <- accdist(d_spatvect = R2025_grid, s_spatvect = f_latr24, roads = roads)
    Dmat_pumps24 <- accdist(d_spatvect = R2025_grid, s_spatvect = f_pumps24, roads = roads)
    Dmat_shower24 <- accdist(d_spatvect = R2025_grid, s_spatvect = f_shower24, roads = roads)
    
    ACC25 <- acc(d_spatvect = R2025_grid, d_attr = "Total25Jan", s_spatvect = f_latr24, s_attr = "nb_Latrine", acolname = "LT_t",distmat = Dmat_latr24)|>
        acc(d_spatvect = _, d_attr = "Total25Jan", s_spatvect = f_pumps24, s_attr = "nb_TW", acolname = "P_t",Dmat_pumps24)|>
        acc(d_spatvect = _, d_attr = "Total25Jan", s_spatvect = f_shower24, s_attr = "nb_WR", acolname = "S_t",Dmat_shower24)|>
        mutate(AVG_t = (LT_t+P_t+S_t)/3)
    
    # free the variable for easing memory pressure
    rm(Dmat_latr24, Dmat_pumps24, Dmat_shower24)
    
    
    cat("Calculating road network distance. Takes some time (Approx. 10 min).\nPlease be patient.\n\n\n")
    # Takes about 20 minutes###
    Dmat_latr22 <- accdist(d_spatvect = R2022_grid, s_spatvect = f_latr22,roads = roads)
    Dmat_pumps22 <- accdist(d_spatvect = R2022_grid, s_spatvect = f_pumps22,roads = roads)
    Dmat_shower22 <- accdist(d_spatvect = R2022_grid, s_spatvect = f_shower22,roads = roads)
    
    
    ACC22 <- acc(d_spatvect = R2022_grid, d_attr = "Total22Feb", s_spatvect = f_latr22, s_attr = "LT", acolname = "LT_t",distmat = Dmat_latr22)|>
        acc(d_spatvect = _, d_attr = "Male22Feb", s_spatvect = f_latr22, s_attr = "LT_Male_sum", acolname = "LT_m",distmat = Dmat_latr22)|>
        acc(d_spatvect = _, d_attr = "Female22Feb", s_spatvect = f_latr22, s_attr = "LT_Female_sum", acolname = "LT_f",distmat = Dmat_latr22)|>
        acc(d_spatvect = _, d_attr = "Total22Feb", s_spatvect = f_pumps22, s_attr = "TW", acolname = "P_t",distmat = Dmat_pumps22)|>
        acc(d_spatvect = _, d_attr = "Total22Feb", s_spatvect = f_shower22, s_attr = "Bathing", acolname = "S_t",distmat = Dmat_shower22)|>
        acc(d_spatvect = _, d_attr = "Male22Feb", s_spatvect = f_shower22, s_attr = "Bathing_M", acolname = "S_m",distmat = Dmat_shower22)|>
        acc(d_spatvect = _, d_attr = "Female22Feb", s_spatvect = f_shower22, s_attr = "Bathing_F", acolname = "S_f",distmat = Dmat_shower22) |>
        mutate(AVG_t = (LT_t+P_t+S_t)/3)

    writeVector(ACC22,"./out/ACC22.gpkg",overwrite = T)
    writeVector(ACC25,"./out/ACC25.gpkg",overwrite = T)
    writeVector(get_ACC_D(acc25 = ACC25, acc22 = ACC22, aoigrid = aoigrid),"./out/ACC_D.gpkg",overwrite = T)
    
    

    
    
    ACC25_survay_testNplot <- function() {
        
        WASH_report <- tryCatch({
            
            get_2024_latr()
            
        },error = function(e){
            
            read_csv("./data/facility/Overview-and-Monitoring-of-WASH-Per-Camp_Round_5_October-31_2024.csv")
            
        })
        
        
        
        ACC25_survey_join <- intersect(ACC25,A1pop) |> as_tibble() |>
            group_by(SMSDCamp) |> summarise(acc25_A1mean = mean(LT_t, na.rm = T) ) |> 
            left_join(select(WASH_report, c("SMSDCamp","people_per_functional_latrine")))
        
        cortest <-  cor.test(ACC25_survey_join$acc25_A1mean, ACC25_survey_join$people_per_functional_latrine, method = "spearman")
        
        utils::capture.output(cortest, file = "./out/FigureS5_cortest.txt")
        
        
        library(ggplot2)
        library(ggtext)
        rho <- as.numeric(cortest$estimate)
        pval <- cortest$p.value
        
        p <- ggplot(ACC25_survey_join,
                    aes(x = acc25_A1mean, y = people_per_functional_latrine)) +
            geom_point(size = 2.2, color = "#1f78b4", alpha = 0.8) +
            geom_smooth(method = "lm", color = "#e31a1c", se = FALSE, linewidth = 1) +
            annotate(
                "text", x = Inf, y = Inf,
                label = sprintf("Spearman ρ = %.2f\np = %.1e ", rho, pval),
                hjust = 1.1, vjust = 1.3, size = 5, family = "sans"
            ) +
            labs(
                #title = "Relationship Between Accessibility and Latrine Pressure",
                x = "Our accessibility scores to latrines (2025)",
                y = "People per functional latrine (2024)"
            ) +
            theme_minimal(base_family = "sans", base_size = 11) +
            theme(
                #plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 9),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_line(linewidth = 0.3, color = "grey80"),
                plot.margin = margin(10, 10, 10, 10)
            )
        
        
        ggsave(
            filename = "./out/FigureS5_validation.pdf",
            plot = p,
            width = 5.91, height = 4.92, units = "in",
            scale = 1,
            dpi = 300,
            device = cairo_pdf,
            bg = "white"
        )
    }
    
    ACC25_survay_testNplot()

    
    # free the variable for easing memory pressure
    rm(Dmat_pumps22)
    
    # Sinario II
    cat("Calculating Scenario II ACC (2/3)\n\n\n")
    f_latr22_S2 <- mutate(f_latr22, LT_Female_sum = LT_all_gen*.75 + LT_Female)
    f_shower22_S2 <- mutate(f_shower22, Bathing_F = Bathing_F*.75 + Bath_gen_u)
    
    acc(d_spatvect = R2022_grid, d_attr = "Total22Feb", s_spatvect = f_latr22_S2, s_attr = "LT", acolname = "LT_t",distmat = Dmat_latr22)|>
        acc(d_spatvect = _, d_attr = "Male22Feb", s_spatvect = f_latr22_S2, s_attr = "LT_Male_sum", acolname = "LT_m",distmat = Dmat_latr22)|>
        acc(d_spatvect = _, d_attr = "Female22Feb", s_spatvect = f_latr22_S2, s_attr = "LT_Female_sum", acolname = "LT_f",distmat = Dmat_latr22)|>
        acc(d_spatvect = _, d_attr = "Total22Feb", s_spatvect = f_shower22_S2, s_attr = "Bathing", acolname = "S_t",distmat = Dmat_shower22)|>
        acc(d_spatvect = _, d_attr = "Male22Feb", s_spatvect = f_shower22_S2, s_attr = "Bathing_M", acolname = "S_m",distmat = Dmat_shower22)|>
        acc(d_spatvect = _, d_attr = "Female22Feb", s_spatvect = f_shower22_S2, s_attr = "Bathing_F", acolname = "S_f",distmat = Dmat_shower22) |>
    writeVector("./out/ACC22_S2.gpkg",overwrite = T)
    
    
    
    # 2SFCA using euclidean distance matrix 
    cat("Calculating ACC using Euclidean distance matrix (3/3)\n\n\n")
    # Takes about less then 1 minute###
    Dmat_latr22 <- accdist(d_spatvect = R2022_grid, s_spatvect = f_latr22)
    Dmat_pumps22 <- accdist(d_spatvect = R2022_grid, s_spatvect = f_pumps22)
    Dmat_shower22 <- accdist(d_spatvect = R2022_grid, s_spatvect = f_shower22)
    

    ############################
    
    ACC22 <- acc(d_spatvect = R2022_grid, d_attr = "Total22Feb", s_spatvect = f_latr22, s_attr = "LT", acolname = "LT_t",distmat = Dmat_latr22)|>
        acc(d_spatvect = _, d_attr = "Male22Feb", s_spatvect = f_latr22, s_attr = "LT_Male_sum", acolname = "LT_m",distmat = Dmat_latr22)|>
        acc(d_spatvect = _, d_attr = "Female22Feb", s_spatvect = f_latr22, s_attr = "LT_Female_sum", acolname = "LT_f",distmat = Dmat_latr22)|>
        acc(d_spatvect = _, d_attr = "Total22Feb", s_spatvect = f_pumps22, s_attr = "TW", acolname = "P_t",distmat = Dmat_pumps22)|>
        acc(d_spatvect = _, d_attr = "Total22Feb", s_spatvect = f_shower22, s_attr = "Bathing", acolname = "S_t",distmat = Dmat_shower22)|>
        acc(d_spatvect = _, d_attr = "Male22Feb", s_spatvect = f_shower22, s_attr = "Bathing_M", acolname = "S_m",distmat = Dmat_shower22)|>
        acc(d_spatvect = _, d_attr = "Female22Feb", s_spatvect = f_shower22, s_attr = "Bathing_F", acolname = "S_f",distmat = Dmat_shower22) |>
        mutate(AVG_t = (LT_t+P_t+S_t)/3)
    
    # free the variable for easing memory pressure
    rm(Dmat_latr22, Dmat_pumps22, Dmat_shower22)
    
    Dmat_latr24 <- accdist(d_spatvect = R2025_grid, s_spatvect = f_latr24)
    Dmat_pumps24 <- accdist(d_spatvect = R2025_grid, s_spatvect = f_pumps24)
    Dmat_shower24 <- accdist(d_spatvect = R2025_grid, s_spatvect = f_shower24)
    
    
    ACC25 <- acc(d_spatvect = R2025_grid, d_attr = "Total25Jan", s_spatvect = f_latr24, s_attr = "nb_Latrine", acolname = "LT_t",distmat = Dmat_latr24)|>
        acc(d_spatvect = _, d_attr = "Total25Jan", s_spatvect = f_pumps24, s_attr = "nb_TW", acolname = "P_t",Dmat_pumps24)|>
        acc(d_spatvect = _, d_attr = "Total25Jan", s_spatvect = f_shower24, s_attr = "nb_WR", acolname = "S_t",Dmat_shower24)|>
        mutate(AVG_t = (LT_t+P_t+S_t)/3)
    
    
    writeVector(ACC22,"./out/ACC22_euclidean.gpkg",overwrite = T)
    writeVector(ACC25,"./out/ACC25_euclidean.gpkg",overwrite = T)
    writeVector(get_ACC_D(acc25 = ACC25, acc22 = ACC22, aoigrid = aoigrid),"./out/ACC_D_euclidean.gpkg",overwrite = T)
    
    cat("\n\n\nAll Done!\n\n\n")
    cat("\n\n\nThe output files of running this script are located at './out/*' \n\n\n")
    
}




if (!interactive()) {
    
    main()
}