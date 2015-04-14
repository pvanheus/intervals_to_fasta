# Taken from Galaxy - included here so we don't need to depend on 
# Galaxy code
import types

def parse_gff_attributes( attr_str ):
    """
    Parses a GFF/GTF attribute string and returns a dictionary of name-value 
    pairs. The general format for a GFF3 attributes string is 
        name1=value1;name2=value2
    The general format for a GTF attribute string is 
        name1 "value1" ; name2 "value2"
    The general format for a GFF attribute string is a single string that
    denotes the interval's group; in this case, method returns a dictionary 
    with a single key-value pair, and key name is 'group'
    """    
    attributes_list = attr_str.split(";")
    attributes = {}
    for name_value_pair in attributes_list:
        # Try splitting by '=' (GFF3) first because spaces are allowed in GFF3
        # attribute; next, try double quotes for GTF.
        pair = name_value_pair.strip().split("=")
        if len( pair ) == 1:
            pair = name_value_pair.strip().split("\"")
        if len( pair ) == 1:
            # Could not split for some reason -- raise exception?
            continue
        if pair == '':
            continue
        name = pair[0].strip()
        if name == '':
            continue
        # Need to strip double quote from values
        value = pair[1].strip(" \"")
        attributes[ name ] = value
        
    if len( attributes ) == 0:
        # Could not split attributes string, so entire string must be 
        # 'group' attribute. This is the case for strictly GFF files.
        attributes['group'] = attr_str
    return attributes

def gff_attributes_to_string(attributes):
    attr_str = ''
    # ensure ID is the first attribute
    if 'ID' in attributes:
        attr_str += 'ID={}; '.format(attributes['ID'])
        del attributes['ID']
    for (key, value) in attributes.iteritems():
        attr_str += '{}={}; '.format(key, value)
    return attr_str.rstrip() # strip off trailing ' '

def gff_string_from_list(fields):
    assert len(fields) == 9, "List to convert to GFF3 must have 9 fields"
    gff_string = ''
    for field in fields[:8]:
         gff_string += str(field) + '\t'
    if type(fields[8]) == types.DictType:
        gff_string += gff_attributes_to_string(fields[8])
    else:
        gff_string += fields[8]
    gff_string += '\n'
    return gff_string
