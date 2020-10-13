import json

def parse_parameter_notebook(notebook_json):

    spike_in_cells = list()
    suffixes = list()
    for cell in notebook_json['cells']:
        if cell['cell_type'] == 'code':
            spike_in_cell = {'cell_type': 'code', 'execution_count': 0, 'metadata': {}, 
                'outputs': [], 'source': cell['source']}
            spike_in_cells.append(spike_in_cell)
        if cell['cell_type'] == 'raw':
            name = cell['source'][0].split()
            assert len(name) == 1
            suffixes.append(name[0])

    if not suffixes or not len(spike_in_cells) == len(suffixes):
        suffixes = list(map(str, range(len(spike_in_cells))))

    return spike_in_cells, suffixes

def extract_parameter_code(file_name, tag):

    with open(file_name) as f:
        notebook_json = json.loads(f.read())

    spike_in_cells, suffixes = parse_parameter_notebook(notebook_json)

    paramter_sets = dict(zip(suffixes, [''.join(cell['source']) for cell in spike_in_cells]))

    return ''.join(paramter_sets[tag])