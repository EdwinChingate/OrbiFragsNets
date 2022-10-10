from IPython.display import HTML, display
import tabulate	
def ShowDF(DF,col=''):
    if col=='':
        col=list(DF.columns)
    display(HTML(tabulate.tabulate(DF[col], headers= col,tablefmt='html')))
    
