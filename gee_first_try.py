# conda activate C:\Users\japolo\Documents\code\condas\gee_env
import geetools
#import folium
import ee

ee.Authenticate()


ee.Initialize()

# Some prelims:
# Outline of lower 48 U.S.
usbound = ee.Image("users/japolo/us_bounds")
usb = usbound.geometry()

#oeel = require('users/OEEL/lib:loadAll')
def clip_img(img):
    return ee.Image(img).clip(usb)

# Set up Daymet collection from 2009 to 2020
# Used just a couple of days to make the test download move faster. Each day takes a couple of minutes to download!
daymetb = ee.ImageCollection("NASA/ORNL/DAYMET_V4").filterDate('2009-06-01', '2009-06-03').map(clip_img)

# Set up projection for reproject of Daymet
proj = ee.Projection('EPSG:4326')

# Wrong
def imreproj(image):
    return image.reproject(crs='EPSG:4326', scale = 1000)

# Reproject Daymet
daymet = daymetb.map(imreproj)

# Vapor pressure from Daymet
dayvp = daymet.select('vp')

# Daylength
daydl = daymet.select('dayl')

# Set a threshold of temps for use in a later calculation
def threshold(image):
    return image.where(image.lt(-5.0), -5)

# Set up tmax and tmin
daytm = daymet.select('tmax','tmin')
# Apply threshold
daytm2 = daytm.map(threshold)

# Use an equals filter to define how the collections match.
filter = ee.Filter.equals(leftField = 'system:index', rightField = 'system:index')

# Create the join.
simpleJoin = ee.Join.simple()

# Applt join
mod1join = ee.ImageCollection(simpleJoin.apply(daytm2, daydl, filter))
mod2join = ee.ImageCollection(simpleJoin.apply(daydl, daytm2, filter))

def enhance(image):
    # Create a collection with 1 image
    temp = ee.ImageCollection(ee.List([image]))
    # Apply join to collection 2
    # Resulting collection will have 1 image with exact same date as img
    join = simpleJoin.apply(mod2join, temp, filter)
    # Get resulting image
    i2 = ee.Image(join.first())
    return image.addBands(i2)

step1 = mod1join.map(enhance)

TC = 4.0
P = 1.5

def temB4sunris(image):
  day_len = image.select('dayl').divide(3600)
  up_time = day_len.multiply(-0.5).add(12)
  set_time = day_len.multiply(0.5).add(12)
  nitelen = ee.Image().expression(
  'daylen * -1 + 24', {
    'daylen' : day_len
  })
  snsetT = ee.Image().expression(
  '(tmin + (tmax - tmin) * sin(pi * (daylen / (daylen + 2 * P))))', {
    'tmin' : image.select('tmin'),
    'tmax' : image.select('tmax'),
    'daylen' : day_len,
    'pi' : 3.14159,
    'P' : 1.5
    })
  snsetT2 = ee.Image().expression(
  '(tmin - snsetT * exp(-1 * nitelen / TC) + (snsetT - tmin) * exp(-1 * (hr + 24 - sunset) / TC)) / (1 - exp(-1 * nitelen / TC))', {
    'tmin' : image.select('tmin'),
    'snsetT' : snsetT,
    'nitelen' : nitelen,
    'tmax' : image.select('tmax'),
    'hr' : 0,
    'sunset' : set_time,
    'daylen' : day_len,
    'pi' : 3.14159,
    'P' : 1.5,
    'TC' : 4.0
    })
  snsetT2 = snsetT2.rename('T_early')
  return snsetT2#.copyProperties(image, ['system:time_start']).copyProperties(image,['system:index'])

mnitetem09 = step1.map(temB4sunris)

mod1join = ee.ImageCollection(simpleJoin.apply(dayvp, mnitetem09, filter))
mod2join = ee.ImageCollection(simpleJoin.apply(mnitetem09, dayvp, filter))

tearlyvp = mod1join.map(enhance)

# l1 = tearlyvp.toList(180)
# t1 = ee.Image(ee.List(l1).get(170))
# print(t1.bandNames().getInfo())
#


# As temp climbs
def tempclimb(image):
   daylt_len = ee.Image(daylen_cal(image))
   up_time = daylt_len.multiply(-0.5).add(12)
   tempup = ee.Image.expression(
    'tmin + (tmax - tmin) * sin(pi * (hr - sunris) / (daylen + 2 * P))', {
      'tmin' : image.select('tmin'),
      'tmax' : image.select('tmax'),
      'pi' : 3.14159,
      'hr' : i,
      'sunris' : up_time,
      'daylen' : daylt_len,
      'P' = 1.5
    )
   return tempup

#before sunset
def tempdown(image):
   daylt_len = ee.Image(daylen_cal(image))
   up_time = daylt_len.multiply(-0.5).add(12)
   tempdrop = ee.Image.expression(
    'tmin + (tmax - tmin) * sin(pi * (hr - sunris) / (daylen + 2 * P))',{
      'tmin' : image.select('tmin'),
      'tmax' : image.select('tmax'),
      'pi' : 3.14159,
      'hr' : i,
      'sunris' : up_time,
      'daylen' : daylt_len,
      'P' = 1.5
    )
   return tempdown

# before midnight
def tempnite = function(image) {
  daylt_len = ee.Image(daylen_cal(image))
  up_time = daylt_len.multiply(-0.5).add(12)
  set_time = daylt_len.multiply(0.5).add(12)
  nitelen = ee.Image().expression(
  'daylen * -1 + 24', {
    'daylen' : daylt_len
  })
  snsetT = ee.Image().expression(
  '(tmin + (tmax - tmin) * sin(pi * (daylen / (daylen + 2 * P))))', {
    'tmin' : image.select('tmin'),
    'tmax' : image.select('tmax'),
    'daylen' : daylt_len,
    'pi' : 3.14159,
    'P' : 1.5
    })
  tempdark = ee.Image.expression(
  '(tmin - snsetT * exp(-1 * nitelen / TC) + (snsetT - tmin) * exp(-1 * (hr - sunset) / TC)) / (1 - exp(-1 * nitelen / TC))', {
    'tmin' : image.select('tmin'),
    'snsetT' : snsetT,
    'nitelen' : nitelen,
    'tmax' : image.select('tmax'),
    'hr' : i,
    'sunset' : set_time,
    'daylen' : daylt_len,
    'pi' : 3.14159,
    'P' : 1.5,
    'TC' : 4.0
    })
   return tempdark

def relhum(image):
    svp1 = ee.Image().expression(
    '611 * 10**(7.5 * temp / (237.7 + temp))',{  #pascals. use 6.11 for kpa
      'temp' : image.select('T_early')
    })
    svp2 = ee.Image().expression(
    'vp * 100 / svp1', {
      'vp' : image.select('vp'),
      'svp1' : svp1
    })
    #svp3 = svp2.where(svp2.gt(100), 100)
    #rel_h = svp3.where(svp3.lt(0), 0)
    rel_h = svp2.rename('rel_hum')
    return ee.Image(rel_h).copyProperties(image, ['system:time_start']).copyProperties(image,['system:index'])

mnittmrh09 = tearlyvp.map(relhum)

mod1join = ee.ImageCollection(simpleJoin.apply(mnittmrh09, tearlyvp, filter))
mod2join = ee.ImageCollection(simpleJoin.apply(tearlyvp, mnittmrh09, filter))

tearlyvp = mod1join.map(enhance)

testcoll = mnittmrh09.filterDate('2009-06-15')

# Set parameters for downloads.
scale = 1000
name = 'temp_00v'
name_pattern = name+'_{system_date}'
## the keywords between curly brackets can be {system_date} for the date of the
## image (formatted using `date_pattern` arg), {id} for the id of the image
## and/or any image property. You can also pass extra keywords using the `extra`
## argument. Also, numeric values can be formatted using a format string (as
## shown in {WRS_PATH:%d} (%d means it will be converted to integer)
date_pattern = 'yMMMdd' # dd: day, MMM: month (JAN), y: year
folder = 'daymets_temp_hr'
data_type = 'uint8'
#region = usb

# ## Export
tasks = geetools.batch.Export.imagecollection.toDrive(
            collection=tearlyvp,
            folder=folder,
            namePattern=name_pattern,
            scale=scale,
            dataType=data_type,
            datePattern=date_pattern,
            verbose=True,
            maxPixels=int(1e13)
        )

### PROBABLY DELETE
lat, lon = 37, -100

my_map = folium.Map(location=[lat,lon], zoom_start=9)
my_map
       
def add_ee_layer(self, ee_image_object, vis_params, name):
    """Adds a method for displaying Earth Engine image tiles to folium map."""
    map_id_dict = ee.Image(ee_image_object).getMapId(vis_params)
    folium.raster_layers.TileLayer(
        tiles=map_id_dict['tile_fetcher'].url_format,
        attr='Map Data &copy; <a href="https://earthengine.google.com/">Google Earth Engine</a>',
        name=name,
        overlay=True,
        control=True
    ).add_to(self)

# Add Earth Engine drawing method to folium.
folium.Map.add_ee_layer = add_ee_layer

vi_p = {
    'min': 1,'max': 17,
    'palette': ['05450a','086a10', '54a708', '78d203', '009900', 'c6b044',
                'dcd159', 'dade48', 'fbff13', 'b6ff05', '27ff87', 'c24f44',
                'a5a5a5', 'ff6d4c', '69fff8', 'f9ffa4', '1c0dff']
}

my_map.add_ee_layer(usbound, vi_p, 'test')

my_map.add_child(folium.LayerControl())

display(my_map)

###dec0112 = ee.Image("users/japolo/time_raster_2012_01_01")
###dec0112 = dec0112.clip(usb)
###
#### Create list of dates for time series
###ndays = ee.Date('2009-06-30').difference('2009-01-01', 'day')
###
###dates = ee.List.sequence(0,ndays,1)
###
####var make_list = function(n) {
####  return Dates_start_10.advance(n,'day')
####}
####dates = dates.map(make_datelist)
###
#### This is critical for building the hour by daily skeleton. It's clunky, but works and is critical.
###def addmoreims(image):
###    return ee.Image(dec0112)
###
#### Soln from developers to create an ImageCollection out of thin air.
#### Actually will create a list, so it's a two-step process
###def CreateCollection(d1):
###  newdate = ee.Date(d1).millis()
###  name = ee.Date(newdate).format('YYYY-MM-dd H:00:00')
###  intCol = ee.ImageCollection(dec0112).map(addmoreims).toList(2)
###  # Change the ImageCollection to an Image
###  # otherwise, the list created will be a list 
###  # of ImageCollections and that can't be 
###  # turned into an ImageCollection
###  return ee.Image(ee.List(intCol).get(0)).set({name: name, "system:time_start": newdate})

# Change list of Images to ...
#var list_of_images = dates.map(CreateCollection)
# ... ImageCollection. Done. For first few thousand hours...
#var ims_l = ee.ImageCollection(list_of_images)
#listoimage = dates.map(CreateCollection)
#datesIC = ee.ImageCollection(listoimage)


#day_len = image.select('dayl').divide(3600)
#up_time = day_len.multiply(-0.5).add(12)
#set_time = day_len.multiply(0.5).add(12)
#nitelen = ee.Image().expression(
#if i < up_time:
#    temtem = temB4sunris(image, i)
#    return temtem
#else if i > up_time & i < 13.5:
#    tem


##### Calculate hours of daylength at each latitude in an Image
##### There is a dayl band in IC!!!!!
####def daylen_cal(image):
####    # Constants
####    pi = 3.141593
####    rad = pi/180
####    # Latitude of pixel
####    imla = ee.Image.pixelLonLat().select('latitude').clip(usb)
####    # Date --- THIS COULD CHANGE FOR ImageCollection vs. Image
####    imgda = ee.Image(image).date()
####    # Day of year
####    doy = imgda.getRelative('day', 'year').add(1)
####    # Sine of latitude
####    latsn = imla.multiply(rad).sin()
####    # Cosine of latitude
####    latcs = imla.multiply(rad).cos()
####    # Max of sine of declination
####    xsindecl = ee.Number(rad).multiply(23.45).sin()
####    # Sine of declination
####    sndl1 = doy.add(10).divide(365).multiply(pi).multiply(2).cos()
####    sindecl = xsindecl.multiply(-1).multiply(sndl1)
####    # Cosine of declination
####    def cosf():
####        return ee.Image().expression(
####          '(1 - im * im)**0.5', {
####          'im' : sindecl
####      })
####    cosdecl = cosf()
####    # Parts for the day length equation
####    p_a = latsn.multiply(sindecl)
####    p_b = latcs.multiply(cosdecl)
####    p_c = p_a.divide(p_b)
####    # Daylength equation
####    def dayeq():
####      return ee.Image().expression(
####        '12 * (1 + (2 / pi) * atan(im / (im*im+1)**0.5))', {
####        'pi' : pi,
####        'im' : p_c
####      })
####    daylen = dayeq()
####    #return ee.Image(imgda)
####    return ee.Image(daylen).copyProperties(image, ['system:time_start']).copyProperties(image,['system:index']).copyProperties(image,['name']).
####
####lite_len = step1.select('tmin').map(daylen_cal)
####
####mod1join = ee.ImageCollection(simpleJoin.apply(step1, lite_len, filter))
####mod2join = ee.ImageCollection(simpleJoin.apply(lite_len, step1, filter))
####
####print('Joined', mod1join, mod2join)
####
####def enhance(image):
####    # Create a collection with 1 image
####    temp = ee.ImageCollection(ee.List([image]))
####    # Apply join to collection 2
####    # Resulting collection will have 1 image with exact same date as img
####    join = simpleJoin.apply(mod2join, temp, filter)
####    # Get resulting image
####    i2 = ee.Image(join.first())
####    return image.addBands(i2)
####
####step2 = mod1join.map(enhance)

#Download.ImageCollection.toDrive = function(collection, folder, options) {
#  var defaults = {
#      scale: 1000,
#      maxPixels: 1e13,
#      type: 'float',
#      region: null,
#      name: null,
#      crs: null,
#      dateFormat: 'yyyy-MM-dd'
#    }
#    
#    var opt = tools.get_options(defaults, options)
#    
#    var n = collection.size().getInfo()
#    
#    var colList = collection.toList(n)
#    
#    if (!opt.name) {
#      var colID = collection.getInfo()['id'] || ""
#      #colID = colID.replace('/','_')
#      colID = helpers.string.formatTask(colID)
#    }
#    
#    for (var i = 0 i < n i++) {
#      var img = ee.Image(colList.get(i))
#      if (opt.name) {
#        if (opt.name.indexOf('{id}') !== -1) {
#          var iid = img.id().getInfo() || i.toString()
#        }
#        if (opt.name.indexOf('{date}') !== -1) {
#          var idate = img.date().format(opt.dateFormat).getInfo()
#        }
#        var id = helpers.string.format(opt.name, {id:iid, date:idate, position:i.toString()})
#        
#        #id = helpers.string._str(id)
#      } else {
#        var id = img.id().getInfo() || colID+'_image_'+i.toString()
#      }
#      
#      var region = opt.region || img.geometry().bounds().getInfo()["coordinates"]
#      
#      var imtype = IMAGE_TYPES(img, opt.type)
#      
#      id = helpers.string.formatTask(id)
#      
#      var params = {
#        image: daym_09_,
#        description: id,
#        folder: temp,
#        fileNamePrefix: daym_09,
#        region: region,
#        scale: opt.scale,
#        maxPixels: opt.maxPixels}
#      
#      if (opt.crs) {
#        params.crs = opt.crs
#      }
#      
#      Export.image.toDrive(params)
#    }
#  }
#
#exports.Download = Download
#
#batch.Download.ImageCollection.toDrive(mnitetem09, 'temp', {image:'daym_09_tem_mnt', crs: 'EPSG:4326'})
