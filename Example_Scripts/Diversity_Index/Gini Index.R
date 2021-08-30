dom_patch_no = read.csv("Results/Patch_Stats/Mean_sum_total/Patchno_matrix_dom_12-08-19.csv")
View(dom_patch_no)
dom_patch_no <- dom_patch_no[,-1]
dom_patch_no <- subset(dom_patch_no, select=-c(BARE,DEAD,TOTAL))
dom_patch_no[is.na(dom_patch_no)] <- 0

library(diverse)
dom_patch_no_mat <- as.matrix(dom_patch_no)
diversity_patch.dom <- diversity(dom_patch_no_mat, type = 'gs', category_row = FALSE)
diversity_patch.dom

## Transposed
# Domiant patches
dom_patch_no_mat.t <- t(dom_patch_no_mat)
t.dom_patch <- diversity(dom_patch_no_mat.t, type = 'gs', category_row = FALSE)
write.csv(t.dom_patch, "Results/Patch_Stats/Mean_sum_total/Gini_Index/dom_patch.csv")

# Subdominant_patches
sub_patch_no = read.csv("Results/Patch_Stats/Mean_sum_total/Patchno_matrix_sub_12-08-19.csv")
View(sub_patch_no)
sub_patch_no <- subset(sub_patch_no, select=-c(quad_no, BARE,DEAD,TOTAL))
sub_patch_no[is.na(sub_patch_no)] <- 0

sub_patch_no_mat <- as.matrix(sub_patch_no)

## Transposed
sub_patch_no_mat.t <- t(sub_patch_no_mat)
t.sub_patch <- diversity(sub_patch_no_mat.t, type = 'gs', category_row = FALSE)
t.sub_patch
write.csv(t.sub_patch, "Results/Patch_Stats/Mean_sum_total/Gini_Index/sub_patch.csv")

# dominant area
dom_area = read.csv("Results/Patch_Stats/Mean_sum_total/Area_spp_dominant_4-08-19.csv")
View(dom_area)
dom_area <- subset(dom_area, select=-c(quad_no,BARE,DEAD,LITTER, ROCK, mean))
dom_area[is.na(dom_area)] <- 0

dom_area_no_mat <- as.matrix(dom_area)

## Transposed
# Domiant area
dom_area_no_mat.t <- t(dom_area_no_mat)
t.dom_area <- diversity(dom_area_no_mat.t, type = 'gs', category_row = FALSE)
t.dom_area
write.csv(t.dom_area, "Results/Patch_Stats/Mean_sum_total/Gini_Index/dom_area.csv")

# subdominant area
sub_area = read.csv("Results/Patch_Stats/Mean_sum_total/Area_spp_sub_4-08-19.csv")
View(sub_area)
sub_area <- subset(sub_area, select=-c(quad_no,BARE,DEAD,LITTER, ROCK, mean))
sub_area[is.na(sub_area)] <- 0

sub_area_no_mat <- as.matrix(sub_area)

## Transposed
# subdomiant area
sub_area_no_mat.t <- t(sub_area_no_mat)
t.sub_area <- diversity(sub_area_no_mat.t, type = 'gs', category_row = FALSE)
t.sub_area
write.csv(t.sub_area, "Results/Patch_Stats/Mean_sum_total/Gini_Index/sub_area.csv")

## Dominant Shape structure
dom_shp= read.csv("Results/Patch_Stats/Mean_sum_total/Shape_dom_matrix_12-08-19.csv")
View(dom_shp)
dom_shp <- subset(dom_shp, select=-c(quad_no,BARE,DEAD,LITTER, ROCK, mean))
dom_shp[is.na(dom_shp)] <- 0

dom_shp_mat <- as.matrix(dom_shp)

## Transposed
dom_shp_mat.t <- t(dom_shp_mat)
t.dom_shp <- diversity(dom_shp_mat.t, type = 'gs', category_row = FALSE)
t.dom_shp
write.csv(t.dom_shp, "Results/Patch_Stats/Mean_sum_total/Gini_Index/dom_shp.csv")

# subdominant area
sub_shp = read.csv("Results/Patch_Stats/Mean_sum_total/Shape_sub_matrix_12-08-19.csv")
View(sub_shp)
sub_shp <-subset(sub_shp, select=-c(quad_no,BARE,DEAD,LITTER, ROCK, mean))
sub_shp[is.na(sub_area)] <- 0

sub_shp_mat <- as.matrix(sub_shp)

## Transposed
sub_shp_mat.t <- t(sub_shp_mat)
t.sub_shp <- diversity(sub_shp_mat.t, type = 'gs', category_row = FALSE)
t.sub_shp
write.csv(t.sub_shp, "Results/Patch_Stats/Mean_sum_total/Gini_Index/sub_shp.csv")
