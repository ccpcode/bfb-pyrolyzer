import json


def read_json(file, comment='//'):
    """
    Read a JSON file that contains comments. Contents of the file are returned
    as a dictionary. Commented lines in the JSON file are ignored.

    Parameters
    ----------
    file : str
        Path to the JSON file containing model parameters
    comment : str
        Format used for comments. Default is '//'.

    Returns
    -------
    dict
        Dictionary of model parameters
    """
    json_str = ''

    with open(file, 'r') as json_file:
        for line in json_file:
            if comment not in line:
                json_str += line

    json_dict = json.loads(json_str)
    return json_dict
