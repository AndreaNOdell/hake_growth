mesh_30 <- make_mesh(hake_sdmTMB_df_complete_sex, xy_cols = c("X", "Y"), cutoff = 30)
plot(mesh_30)

mesh30 = sdmTMB( 
  data = hake_sdmTMB_df_complete_sex,
  formula = weight ~ s(new_age) + s(cohort) + catch_month + sex_description,
  mesh = mesh_30, 
  time = "catch_year",
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "ar1",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016, 2018, 2020))

mesh_50 <- make_mesh(hake_sdmTMB_df_complete_sex, xy_cols = c("X", "Y"), cutoff = 50)
mesh50 = sdmTMB( 
  data = hake_sdmTMB_df_complete_sex,
  formula = weight ~ s(new_age) + s(cohort) + catch_month + sex_description,
  mesh = mesh_50, 
  time = "catch_year",
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "ar1",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016, 2018, 2020))

mesh_20 <- make_mesh(hake_sdmTMB_df_complete_sex, xy_cols = c("X", "Y"), cutoff = 20)
mesh20 = sdmTMB( 
  data = hake_sdmTMB_df_complete_sex,
  formula = weight ~ s(new_age) + s(cohort) + catch_month + sex_description,
  mesh = mesh_20, 
  time = "catch_year",
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "ar1",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016, 2018, 2020))
