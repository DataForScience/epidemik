from collections import defaultdict

class NotInitialized(Exception):
    pass

epi_colors = defaultdict(lambda :'#f39019')
epi_colors['S'] = '#51a7f9'
epi_colors['E'] = '#f9e351'
epi_colors['I'] =  '#cf51f9'
epi_colors['R'] = '#70bf41'
epi_colors['D'] = '#8b8b8b'