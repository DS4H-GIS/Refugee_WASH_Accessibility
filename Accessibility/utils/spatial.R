library(sf)
library(dplyr)
library(terra)
library(tidyterra)

vectNreproject <- function(dir,terra_crs,bound = NA) {
    
    
    ret <- vect(dir) |> project(terra_crs)
    
    return(ret)
    
    # if(!is.na(bound)) return(ret)
    # 
    # return(ret[!is.na(ret$area)])
    
}


vect_make_grid <- function(bound, cellsize = c(50,50), clip = F){

    copy_bound <- bound
    if (inherits(copy_bound, "SpatVector")) copy_bound <- vect(ext(copy_bound),crs(bound))

    if (!any(c("sf", "sfc") %in% attributes(copy_bound)$class)) copy_bound <- st_as_sf(copy_bound)

    aoigrid <- st_make_grid(copy_bound, cellsize = cellsize) |>
        st_as_sf() |> mutate(g_index = row_number()) |> vect()

    if (clip) return(aoigrid[bound])

    return(aoigrid)

}


# for future ideas
# vect_make_grid <- function(bound, cellsize = c(50, 50), clip = FALSE) {
# 
#     # bound를 SpatVector로 변환 (벡터 객체가 아닌 경우)
#     if (!inherits(bound, "SpatVector")) {
#         bound <- vect(bound)
#     }
# 
#     # 바운드의 extent 얻기
#     ext <- ext(bound)
# 
#     # extent와 cellsize를 이용해 raster 생성
#     r <- rast(ext,crs = crs(bound), resolution = cellsize)
# 
#     # 생성된 raster를 폴리곤으로 변환
#     aoigrid <- as.polygons(r)
# 
#     # 그리드의 인덱스를 추가
#     aoigrid$g_index <- 1:nrow(aoigrid)
# 
#     # clip 인자에 따라 결과 반환
#     if (clip) {
#         return(intersect(aoigrid, bound))  # intersect로 클리핑
#     }
# 
#     return(aoigrid)
# }






# st_intersection(R2022,aoigrid) |>
#     #as.data.table()|>
#     group_by(g_index) |> summarise(area = sum(area))
#     #_[,]


# terra::intersect is 1.7 timesfaster
# tic();st_intersection(R2022,aoigrid);toc()
# tic();intersect(v1,v2);toc()

weighted_sum <- function(source_vect,target_vect,aux_vect = NULL, aux_attr = NULL) {
    
    
    
    if(is.null(aux_vect)){
        
        vect_intersect <- intersect(target_vect,source_vect)
        vect_intersect |> 
            mutate(area = expanse(vect_intersect)) |> 
            select(area) |>
            terra::intersect(target_vect) |>
            as_tibble() |> group_by(g_index) |>
            summarise(area = sum(area,na.rm = T)) |> 
            merge(x=target_vect,y=_,
                  all.x=TRUE
            ) -> ret
        
        return(ret[!is.na(ret$area)])
        
    } else {
        
        
        first_intersect <- intersect(source_vect, aux_vect) |> select(-area)
        first_intersect <- mutate(first_intersect,area = expanse(first_intersect))
        
        first_intersect |> 
            as_tibble() |>
            group_by(CampLabel) |>
            summarise(sum_area = sum(area,na.rm = T)) |>
            merge(x=first_intersect,y=_, all.x=TRUE) |>
            mutate(prop = area/sum_area) |> 
            mutate(across(all_of(aux_attr), ~ .x * prop)) -> first_intersect
        
        
        second_intersect <- intersect(target_vect, first_intersect) |> select(-c(sum_area,prop))
        second_intersect <- mutate(second_intersect, area_split = expanse(second_intersect))
        
        
        second_intersect |> 
            select(all_of(aux_attr), area, area_split) |>
            mutate(prop = area_split/area) |> 
            mutate(across(all_of(aux_attr), ~ .x * prop)) |>
            terra::intersect(target_vect) |>
            as_tibble() |>
            group_by(g_index) |>
            summarise(across(all_of(aux_attr), ~ sum(.x, na.rm = TRUE)))|> 
            merge(x=target_vect,y=_,
                  all.x=TRUE) -> ret
        
        
        return(drop_na(ret))
        
    }
    

    
    
    
} 
