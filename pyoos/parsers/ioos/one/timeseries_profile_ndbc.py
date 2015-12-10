from dateutil import parser
from collections import OrderedDict, defaultdict
from bisect import bisect_left
from urllib import unquote
from copy import copy

from owslib.crs import Crs
from owslib.namespaces import Namespaces
from owslib.util import nspath_eval

from shapely.geometry import Point as sPoint

from paegan.cdm.dsg.member import Member
from paegan.cdm.dsg.features.station_profile import StationProfile
from paegan.cdm.dsg.features.base.point import Point
from paegan.cdm.dsg.features.base.profile import Profile
from paegan.cdm.dsg.collections.base.profile_collection import ProfileCollection
from paegan.cdm.dsg.collections.station_collection import StationCollection

from owslib.swe.common import Time, DataChoice, DataRecord


def get_namespaces():
    ns = Namespaces()
    return ns.get_namespaces(["swe20"])
namespaces = get_namespaces()

def nspv(path):
    return nspath_eval(path, namespaces)

class ProfileCache(object):
    """
    Helper class to accumulate observations and transform them into ProfileCollections representing TimeseriesProfiles.

    Used internally.
    """
    def __init__(self):
        self._cache = defaultdict(OrderedDict)

    def add_obs(self, sensor, t, obs_point, obs_members):
        """
        """
        profile = self._get_profile(sensor, t)
        point   = self._get_point(profile, obs_point)
        point.members.extend(obs_members)

    def get_collections(self):
        return {k[2] : ProfileCollection(elements=pd.values()) for k, pd in self._cache.iteritems()}

    def _get_profile(self, sensor, t):
        sens_loc_tuple = (sensor['location']['point'].x, sensor['location']['point'].y, sensor['station'])
        profile_od     = self._cache[sens_loc_tuple]

        if t not in profile_od:
            profile          = Profile()
            profile.location = sPoint(*sens_loc_tuple[0:2])
            profile.time     = t

            # @TODO this is a hack until we can figure out how to assoc stations properly
            profile.station  = sensor['station']

            profile_od[t]    = profile
            return profile

        return profile_od[t]

    def _get_point(self, profile, point):
        """
        Finds the given point in the profile, or adds it in sorted z order.
        """
        cur_points_z = [p.location.z for p in profile.elements]
        try:
            cur_idx = cur_points_z.index(point.z)
            return profile.elements[cur_idx]
        except ValueError:
            new_idx            = bisect_left(cur_points_z, point.z)
            new_point          = Point()
            new_point.location = sPoint(point)
            new_point.time     = profile.time
            profile.elements.insert(new_idx, new_point)
            return new_point


class TimeSeriesProfileNdbc(object):
    def __init__(self, element):
      # starting iwht data steam element
        elemType = element.find(nspv("swe20:elementType"))
        if elemType is None:
            return
        
        dataRec = elemType.find(nspv("swe20:DataRecord"))
        if dataRec is None:
            return
        
        record = DataRecord(dataRec)

        columns = []
        data_array_columns = {}
        s = StationProfile()
        sensors = {}
        varUnits={}
        varStandardNames={}
        constant_time = None
        for idx, field in enumerate(record.field):
            
            
            #if field.content.value == None:
            #    value_list_order.append(field)
            if field.name == "stationID":
                s.uid   = field.content.value
                s.name = s.uid
            elif field.name == "location":
                vector  = record.get_by_name("location").content
                srss = vector.referenceFrame
                hsrs = None
                try:
                    hsrs = Crs(srss[0])
                except ValueError:
                    pass
                s.set_property("horizontal_srs", hsrs)
                s.set_property("vertical_srs",   None)
                s.set_property("localFrame",    None)

                lat = vector.get_by_name("latitude").content.value
                lon = vector.get_by_name("longitude").content.value
                height = vector.get_by_name("height")
                z = None
                if height is not None:
                    z   = height.content.value
                loc = [lon, lat]
                if z:
                    loc.append(z)
                s.location = sPoint(*loc) 
            elif field.name == "sensorID":
                if field.content.value == None:
                    columns.append(field)
            elif field.name == "time" or field.name == 'date_time':
                if field.content.value == None:
                    columns.append(field)
                else:
                    constant_time = field.content.value    
            elif field.name == "altitude":
                height      = None
                altitude_quantity = record.get_by_name("altitude").content
                if altitude_quantity.value == None:
                    columns.append(field)
                elif altitude_quantity.referenceFrame == "#%s_frame" % s.name:
                    # Uses the station as reference frame
                    if z:
                        height      = z + altitude_quantity.value
                        horizontal_srs  = s.get_property("horizontal_srs")
                        vertical_srs    = s.get_property("vertical_srs")
                else:
                    height      = altitude_quantity.value
                    horizontal_srs  = None
                    vertical_srs    = None
                    if hasattr(altitude_quantity, 'referenceFrame'):
                        srss            = altitude_quantity.referenceFrame
                        try:
                            horizontal_srs = Crs(srss[0])
                        except ValueError:
                            pass
                        try:
                            vertical_srs = Crs(srss[-1].replace("2=http:", "http:"))
                        except ValueError:
                            pass
                        
                        loc = [s.location.x, s.location.y]
                        if height:
                            loc.append(height)

                       

            else:
                try:
                    if field.content.value == None:
                        columns.append(field)
                    varStandardNames[field.name] = field.definition
                    varUnits[field.name] = field.uom
                except:
                
                    # check to see if dealing with a data array
                    if isinstance(field.content, DataArray):
                        element_count_field = field.content.elementCountParent
                        count_field = field.content.elementCount
                        if count_field.text == None:
                            columns.append(element_count_field)
                        element_count_name = element_count_field.attrib['name']
                        element_type = field.content.elementType
                        data_array_columns[element_count_field] = []
                        for field in element_type.content.field:
                            if field.content.value == None:
                                data_array_columns[element_count_field].append(field)


        # Start building the column structure
        data_array =  element.find(nspv("swe20:values"))
        encoding = element.find(nspv("swe20:encoding"))
        t_encoding = encoding.find(nspv("swe20:TextEncoding"))
        extracted_data = {}
        for column in columns:
            extracted_data[column] = []

        location            = sPoint(*loc)
        decimalSeparator    = t_encoding.attrib['decimalSeparator']
        tokenSeparator      = t_encoding.attrib['tokenSeparator']
        blockSeparator      = t_encoding.attrib['blockSeparator']
       
        blockSeparator = unquote(blockSeparator)
        tokenSeparator = unquote(tokenSeparator)
        pt_list = []
        data_values = data_array.text
        self.raw_data = copy(data_values)
        all_depths = set()
        for row in filter(lambda x: x != "", data_values.split("\n")):

            pt = None
            members     = []
            values      = row.split(tokenSeparator)
            sensor_key  = None
            sensor_name = None
            value = None
            var_name = None
            time_index = None
            skip_entries = 0
            i = 0
            curr_loc = location
            for x in columns:
        
                if data_array_columns.has_key(x):
                    data_array_count = int(values[i])
                    num_entries = len(data_array_columns[x])
                    skip_entries = data_array_count * num_entries
                    i += skip_entries
                elif x.name == 'time':
                    time_index = i
                elif x.name == 'altitude':
                    curr_loc = sPoint(*(location.x, location.y, float(values[i])))
                    all_depths.add(float(values[i]))
                elif x.name == 'sensorID':
                    sensor_name =  values[i]
                elif data_array_columns.has_key(x):
                    data_array_count = values[i]
                    num_entries = len(data_array_columns[x])
                    skip_entries = data_array_count * num_entries
                else:
                    value = values[i]
                    var_name = x.name
                    pt      = Point()
                    if constant_time is not None:
                        pt.time = constant_time
                    else:
                        pt.time = parser.parse(values[time_index])

                    m = Member( units=varUnits[var_name],
                                        name=var_name,
                                        standard=varStandardNames[var_name],
                                        value=value)
                    members.append(m)
            
                i += 1
            pt.members = members
            pt.location = curr_loc
            pt_list.append(pt)
            
            #sensors[sensor_key]['values'].append(pt)

#         for k, v in stations.iteritems():
#             for sk, sv in sensors.iteritems():
#                 # Match on station uid
#                 if sv['station'] == k:
#                     v.elements = self._merge_points(v.elements or [], sv['values'])
#  
#         if len(stations) > 1:
#             self.feature = StationCollection(elements=stations)
#         elif len(stations) == 1:
        profileCollection = ProfileCollection(elements=pt_list)
        profileCollection.set_depth_range(all_depths)
        self.feature = profileCollection
        #else:
         #   print "No stations found!"

    def _parse_sensor_orientation(self, ori_el):
        # 'srs':Crs(),    # @TODO (OWSLib cannot parse this Crs yet)
        orientation = {}
        for coord in ori_el.content.coordinate:
            orientation[coord.axisID] = {
                'name' : coord.name,
                'value' : coord.value,
                'axis' : coord.axisID,
                'uom' : coord.uom
            }

        return orientation

    def _parse_location(self, loc_el, station_point):
        vector         = loc_el.content

        srss           = vector.referenceFrame.split("&amp;")

        hsrs = None
        try:
            hsrs = Crs(srss[0])
        except ValueError:
            pass

        vsrs = None
        try:
            vsrs = Crs(srss[-1].replace("2=http:", "http:"))
        except ValueError:
            pass

        local_frame    = vector.localFrame
        lat            = vector.get_by_name("latitude").content.value
        lon            = vector.get_by_name("longitude").content.value
        z              = vector.get_by_name("height").content.value

        loc = [lon, lat]
        if z:
            loc.append(z)
        else:
            loc.append(station_point.z)

        location = {'horizontal_srs' : hsrs,
                    'vertical_srs'   : vsrs,
                    'localFrame'     : local_frame,
                    'point'          : sPoint(*loc)}

        return location

    def _parse_profile_bins(self, profbins_el):
        data_array = profbins_el.content
        count      = int(data_array.elementCount[0].text)
        data       = self._parse_data_array(data_array)

        bin_center_quantity = data_array.elementType.content.get_by_name('binCenter')
        bin_center = {'referenceFrame' : bin_center_quantity.content.referenceFrame,
                      'axisID'         : bin_center_quantity.content.axisID,
                      'uom'            : bin_center_quantity.content.uom,
                      'values'         : data[0]}

        bin_edges_quantityrange = data_array.elementType.content.get_by_name('binEdges')
        bin_edges = {'referenceFrame'  : bin_edges_quantityrange.content.referenceFrame,
                     'axisID'          : bin_edges_quantityrange.content.axisID,
                     'uom'             : bin_edges_quantityrange.content.uom,
                     'values'          : data[1]}

        profile_bins = {'bin_center': bin_center,
                        'bin_edges': bin_edges}

        return profile_bins

    def _parse_profile_heights(self, profheights_el):
        data_array = profheights_el.content
        count      = int(data_array.elementCount[0].text)
        data       = self._parse_data_array(data_array)

        height_el  = data_array.elementType.get_by_name('height')
        profile_definition = {'referenceFrame': height_el.content.referenceFrame,
                              'axisID'        : height_el.content.axisID,
                              'uom'           : height_el.content.uom,
                              'values'        : data[0]}

        return profile_definition

    def _parse_data_array(self, data_array):
        """
        Parses a general DataArray.
        """
        decimalSeparator    = data_array.encoding.decimalSeparator
        tokenSeparator      = data_array.encoding.tokenSeparator
        blockSeparator      = data_array.encoding.blockSeparator
        collapseWhiteSpaces = data_array.encoding.collapseWhiteSpaces

        data_values = data_array.values
        lines = filter(lambda x: x != '', data_values.split(blockSeparator))

        ret_val = []

        for row in lines:
            values = row.split(tokenSeparator)
            ret_val.append([float(v) if ' ' not in v.strip() else [float(vv) for vv in v.split()] for v in values])

        # transpose into columns
        return [list(x) for x in zip(*ret_val)]

    def _parse_sensor_data(self, obs_el, sensor_info):
        """
        Returns ProfileCollection
        """
        data_array = obs_el.content

        # get column defs
        data_record = data_array.elementType.content
        columns = []
        for f in data_record.field:
            columns.append(f)

        # get more information on sensor cols
        sensor_cols = defaultdict(list)
        sensor_vals = defaultdict(list)

        sensor_rec = data_record.get_by_name('sensor')
        for sendata in sensor_rec.content.item:
            if sendata.content is None:
                continue

            for f in sendata.content.field:
                sensor_cols[sendata.name].append(f)

        # @TODO deduplicate
        decimalSeparator    = data_array.encoding.decimalSeparator
        tokenSeparator      = data_array.encoding.tokenSeparator
        blockSeparator      = data_array.encoding.blockSeparator
        collapseWhiteSpaces = data_array.encoding.collapseWhiteSpaces

        data_values = data_array.values
        lines = filter(lambda x: x != '', data_values.split(blockSeparator))

        # profile cacher!
        profile_cache = ProfileCache()

        for row in lines:
            values = row.split(tokenSeparator)

            skey     = None
            i        = 0
            cur_time = None
            cur_qual = None

            for c in columns:

                if isinstance(c.content, Time) and c.content.definition == "http://www.opengis.net/def/property/OGC/0/SamplingTime":
                    cur_time = parser.parse(values[i])
                    i += 1

                    if len(c.quality):
                        # @TODO: do some quality constraint checks
                        i += len(c.quality)
                        # for qua in c.quality:

                elif isinstance(c.content, DataChoice) and c.name == "sensor":
                    sensor_key = values[i]
                    i += 1

                    sensor_dr = c.content.get_by_name(sensor_key).content
                    sensor_info_ = sensor_info[sensor_key]
                    parsed, nc = self._parse_sensor_record(sensor_dr, sensor_info_, values[i:])

                    # turn these into Points/Members
                    for rec in parsed:
                        # calc a Z value from rec/sensor and build point
                        point, members = self._build_obs_point(sensor_info_, rec)

                        # add to profile
                        profile_cache.add_obs(sensor_info_, cur_time, point, members)

                    i += nc

        return profile_cache.get_collections()

    def _parse_sensor_record(self, sensor_data_rec, sensor_info, rem_values):
        """
        Parses values via sensor data record passed in.
        Returns parsed values AND how many items it consumed out of rem_values.
        """
        val_idx = 0

        # @TODO seems there is only a single field in each of these
        assert len(sensor_data_rec.field) == 1
        sensor_data_array = sensor_data_rec.field[0].content

        # there is probably not going to be a count in the def, it'll be in the data
        count = None
        count_text = sensor_data_array.elementCount.text
        if count_text:
            count = int(count_text.strip())

        if not count:
            count = int(rem_values[val_idx])
            val_idx += 1

        parsed = []

        for recnum in xrange(count):
            cur  = []

            for f in sensor_data_array.elementType.field:
                cur_val = rem_values[val_idx]
                val_idx += 1

                m = Member(name=f.name,
                           standard=f.content.definition)

                if hasattr(f.content, 'uom'):
                    m['units'] = f.content.uom

                try:
                    m['value'] = float(cur_val)
                except ValueError:
                    m['value'] = cur_val

                if len(f.quality):
                    m['quality'] = []
                    for qual in f.quality:
                        cur_qual = rem_values[val_idx]
                        val_idx += 1

                        # @TODO check this against constraints
                        m['quality'].append(cur_qual)

                cur.append(m)

            parsed.append(cur)

        return parsed, val_idx

    def _build_obs_point(self, sensor_info, obs_recs):
        """
        Pulls bin/profile height info out and calcs a z.

        TODO: currently extremely naive

        Returns a 2-tuple: point, remaining obs_recs
        """
        cur_point = sensor_info['location']['point']

        keys = [m['name'] for m in obs_recs]
        if 'binIndex' in keys:
            zidx      = keys.index('binIndex')
            bin_index = int(obs_recs[zidx]['value'])
            z         = sensor_info['profile_heights']['values'][bin_index]

            point     = sPoint(cur_point.x, cur_point.y, cur_point.z + z)

        elif 'profileIndex' in keys:
            zidx      = keys.index('profileIndex')
            bin_index = int(obs_recs[zidx]['value'])

            # @TODO take into account orientation, may change x/y/z
            # @TODO bin edges?
            z         = sensor_info['profile_bins']['bin_center']['values'][bin_index]

            point     = sPoint(cur_point.x, cur_point.y, cur_point.z + z)

        else:
            raise ValueError("no binIndex or profileIndex in Member: %s", keys)

        # remove z related Member
        obs_recs = obs_recs[:]
        obs_recs.pop(zidx)

        return point, obs_recs
