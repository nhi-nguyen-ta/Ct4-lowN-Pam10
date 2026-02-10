# def hello() -> str:
#     return "Hello from genome-data-analysis!"

# print(f"Module {__name__} loaded.")


# The logger could be set in here...

# The import below allows us to import simgle functions from single modules
# for example 'from genome-data-analysis.systools import makedir'
# Additionally this also allows us to use the package name and submosules when calling functions
# For example: import 'genome_data_analysis as gda' and then call gda.systools.makedir()
from . import systools
from . import annotation

# This defined the submodule that should be imported when doing 'from genome_data_analysis import *'
# __all__ = ["systools"]
