################################################################################
#                  Animated hexbin map of nightlight in Europe
#                  Milos Popovic
#                  2021/12/29
################################################################################

if(!require("sf")) install.packages("sf")
if(!require("tidyverse")) install.packages("tidyverse")  
if(!require("grid")) install.packages("grid")  
if(!require("gridExtra")) install.packages("gridExtra") 
if(!require("raster")) install.packages("raster")  
if(!require("exactextractr")) install.packages("exactextractr")  
if(!require("gganimate")) install.packages("gganimate") 
if(!require("transformr")) install.packages("transformr") 
if(!require("tidyr")) install.packages("tidyr")
if(!require("giscoR")) install.packages("giscoR") 
if(!require("classInt")) install.packages("classInt")  
if(!require("gifski")) install.packages("gifski")
if(!require("rfigshare")) install.packages("rfigshare")
if(!require("jsonlite")) install.packages("jsonlite")
      
library(tidyverse, quietly=T) #data wrangling
library(tidyr, quietly=T) #data wrangling
library(sf, quietly=T) #geospatial analysis
library(grid, quietly=T) #make grid
library(gridExtra, quietly=T) #make grid
library(raster, quietly=T) # import raster files
library(exactextractr, quietly=T) # zonal statistics
library(gganimate, quietly=T) # animation
library(giscoR, quietly=T) #shapefile of Europe
library(classInt, quietly=T) #bins
library(gifski, quietly=T) #renderer

windowsFonts(georg = windowsFont('Georgia')) #font

# 1. EUROPE SHAPEFILE
europe <- giscoR::gisco_get_countries(
  year = "2016",
  epsg = "4326",
  resolution = "10",
  region = "Europe"
)

# define longlat projection
crsLONGLAT <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# create a bounding box of Europe based on longlat projection
bb <- st_sfc(
  st_polygon(list(cbind(
    c(-23.6, 47.3, 47.3, -23.6, -23.6), # x-coordinates (longitudes)
    c(30.5, 30.5, 71.05, 71.05, 30.5)     # y-coordinates (latitudes)
    ))), crs = crsLONGLAT)

# crop Europe by bounding box
eur <- st_crop(europe, st_bbox(bb))

# transform the shp into Lambert projection EPSG 3575 with meter units
eur2 <- eur %>% st_transform(3575)

# Make grid 25 km
gr <- st_make_grid(eur2, cellsize = 25000, what = "polygons", square = F) %>%
  st_intersection(eur2) %>%
  st_sf() %>%
  mutate(id = row_number())

gg <- gr %>% filter(
  st_geometry_type(.)
  %in% c("POLYGON", "MULTIPOLYGON")) %>% 
  st_cast("MULTIPOLYGON") #transform to multipolygons

g <- st_transform(gg, crs = crsLONGLAT) #transform back to longlat

# 2. NIGHTLIGHT DATA
# list of all nightlight raster files from FigShare
x = rfigshare::fs_details("9828827", mine = F, session = NULL)
urls = jsonlite::fromJSON(jsonlite::toJSON(x$files), flatten = T)

# download harmonized nightlight data for 1992-2020
urls$download <- paste0(urls$download_url, "/", urls$name) 
urls <- urls[-28,] # get rid of the 2013 duplicate
url <- urls[,5]

for (u in url){
  download.file(u, basename(u), mode="wb")
}

# enlist raster files, merge them into a single file and re-project
rastfiles <- list.files(path = getwd(), 
    pattern = ".tif$",
    all.files=T, 
    full.names=F)

dd <- list()
# compute average nighlight DN values for every hexbin
for(i in rastfiles) dd[[i]] <- {
        allrasters <- raster(i)
        d  <- exact_extract(allrasters, g, "mean")
        d
}

dd <- do.call(rbind, dd) # pull them all together
a <- as.data.frame(t(dd)) #convert to data.frame and transpose
a$id <- 1:max(nrow(a)) # id for joining with sf object
names(a) #check column names
ee <- pivot_longer(a, cols=1:29, names_to = "year", values_to = "value") #wide to long format
l <- ee %>% # extract year from raster name
     mutate_at("year", str_replace_all, c("Harmonized_DN_NTL_"="", "_simVIIRS"="", "_calDMSP"="", ".tif"=""))
head(l) 

f <- left_join(l, gg, by="id") %>% # join with transformed sf object
     st_as_sf() %>%
     filter(!st_is_empty(.)) #remove null features

f$year <- as.numeric(as.character(f$year))

# bins
brk <- round(classIntervals(f$value, 
              n = 6, 
              style = 'fisher')$brks, 0)
vmax <- max(f$value, na.rm=T)
vmin <- min(f$value, na.rm=T)
# define the color palette
cols =c("#182833", "#1f4762", "#FFD966", "#ffc619")
newcol <- colorRampPalette(cols)
ncols <- 7
cols2 <- newcol(ncols)

# animate the map
pp <- ggplot(f) +
  geom_sf(mapping = aes(fill = value), color = NA, size=0) +
  geom_sf(data=eur2, fill=NA, color="white", size=0.01) +
  scale_fill_gradientn(name="DN values",
                       colours=cols2,
                       breaks=brk,
                       labels=brk,
                       limits=c(vmin,vmax)) +
  guides(fill=guide_legend(
            direction = "horizontal",
            keyheight = unit(0.25, units = "mm"),
            keywidth = unit(4, units = "mm"),
            title.position = 'top',
            title.hjust = 0.5,
            label.hjust = .5,
            nrow = 1,
            byrow = T,
            reverse = F,
            label.position = "bottom"
          )
    ) +
  theme_minimal() +
  theme(text = element_text(family = "georg", color = "#22211d"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = c(.5, .005),
    legend.text = element_text(size=3, color="white"),
    legend.title = element_text(size=3, color="white"),
    panel.grid.major = element_line(color = "grey20", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face="bold", size=5, color="white", hjust=.5),
    plot.caption = element_text(size=2, color="white", hjust=.5, vjust=-10),
    plot.subtitle = element_text(size=8, color="#ffc619", hjust=.5),
    plot.margin = unit(c(t=0, r=0, b=0, l=0),"lines"), #added these narrower margins to enlarge maps
    plot.background = element_rect(fill = "grey20", color = NA), 
    panel.background = element_rect(fill = "grey20", color = NA), 
    legend.background = element_rect(fill = "grey20", color = NA),
    panel.border = element_blank())+
  labs(x = "", 
         y = NULL, 
         title = "Average nightlight values (1992-2020)", 
         subtitle = 'Year: {as.integer(frame_time)}', 
         caption = "©2021 Milos Popovic (https://milospopovic.net)\n Data: Li et al. (2020). A harmonized global nighttime light dataset 1992–2018. Scientific data, 7(1), 1-9.")

pp1 = pp + 
  transition_time(year) + 
  ease_aes('linear') +
  enter_fade() +
  exit_fade()

p_anim <- animate(pp1,
  nframes = 140,
  duration = 25,
  start_pause = 5,
  end_pause = 30, 
  height = 2140, 
  width = 2830,
  res = 600)

anim_save("nightlight_1992_2020b.gif", p_anim)