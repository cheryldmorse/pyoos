from dateutil import parser
from copy import copy

from owslib.crs import Crs
from owslib.util import nspath_eval
from owslib.namespaces import Namespaces
from shapely.geometry import Point as sPoint
from urllib import unquote
from paegan.cdm.dsg.member import Member
from paegan.cdm.dsg.features.station import Station
from owslib.swe.common import DataArray

from paegan.cdm.dsg.collections.base.point_collection import PointCollection
from paegan.cdm.dsg.features.base.point import Point

from owslib.swe.common import Time, DataChoice, DataRecord, AbstractSimpleComponent


def get_namespaces():
    ns = Namespaces()
    return ns.get_namespaces(["swe20"])
namespaces = get_namespaces()


def nspv(path):
    return nspath_eval(path, namespaces)


class TimeSeries_ndbc(object):
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
        s = Station()
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
                location_quantity = record.get_by_name("altitude").content
                
                altitude_quantity = record.get_by_name("altitude").content
                s.set_property("altitude_units",altitude_quantity.uom)

                if location_quantity.referenceFrame == "#%s_frame" % s.name:
                    # Uses the station as reference frame
                    if location_quantity.value and z:
                        height      = z + location_quantity.value
                        horizontal_srs  = s.get_property("horizontal_srs")
                        vertical_srs    = s.get_property("vertical_srs")
                else:
                    # Uses its own height
                    if location_quantity.value:
                        height      = location_quantity.value
                    horizontal_srs  = None
                    vertical_srs    = None
                    if hasattr(location_quantity, 'referenceFrame'):
                        srss            = location_quantity.referenceFrame
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

                        location            = sPoint(*loc)

            else:
                try:
                    if field.content.value == None:
                        columns.append(field)
                    varStandardNames[field.name] = field.definition
                    varUnits[field.name] = field.uom
                except:
                
                    # check to see if dealing with a data array
                    if isinstance(field.content, DataArray):
                        if hasattr(field.content, 'elementCountParent'):
                            element_count_field = field.content.elementCountParent
                        else:
                            element_count_field = field.content.elementCount

                        count_field = field.content.elementCount
                        if count_field.text == None:
                            columns.append(element_count_field)
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

        decimalSeparator    = t_encoding.attrib['decimalSeparator']
        tokenSeparator      = t_encoding.attrib['tokenSeparator']
        blockSeparator      = t_encoding.attrib['blockSeparator']
       
        blockSeparator = unquote(blockSeparator)
        tokenSeparator = unquote(tokenSeparator)
        pt_list = []
        data_values = data_array.text
        self.raw_data = copy(data_values)
        
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
            for x in columns:
        
                if data_array_columns.has_key(x):
                    data_array_count = int(values[i])
                    num_entries = len(data_array_columns[x])
                    skip_entries = data_array_count * num_entries
                    i += skip_entries
                elif x.name == 'time':
                    time_index = i

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
                                        sensor=sensor_name,
                                        value=value)
                    members.append(m)
            
                i += 1
            pt.members = members
            pt.location = location
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
        self.feature = PointCollection(elements=pt_list)
        self.depth_units = s.get_property("altitude_units")

        #else:
         #   print "No stations found!"

    def _merge_points(self, pc1, pc2):
        """
        Merges points based on time/location.

        @TODO: move to paegan, SO SLOW
        """
        res = pc1[:]

        for p in pc2:
            for sp in res:
                if sp.time == p.time and (sp.location is None or (sp.location.equals(p.location))):
                    sp.members.extend(p.members)
                    break
            else:
                res.append(p)

        return res

