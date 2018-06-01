import math
import time
from clawpack.geoclaw.data import LAT2METER

"""Unit Converter

@param param1: category (l [length], p [pressure_rigidity], s [speed], m [moment], t [temperature])
@param param2: unit from
@param param3: unit to
@param param4: value
@return: converted value

"""

def convertUnits(cat, unit1, unit2, num1):  #add docstring
#def convertUnits():
    Length_conv_factor= {
        'm': 1.0,
        'cm':1e-2,
        'km':1e3,
        'miles': 1.60934e3,
        'nm':1852.0,
        'lat-long': LAT2METER,
    }

    Pressure_Rigidity_conv_factor= {
        'Pa': 1.0,
        'kPa':1e3,
        'mPa':1e6,
        'GPa': 1.e9,
        'mbar':1e2,
        'dyne/cm^2': 0.1,
        'dyne/m^2': 1.e-5,
    }

    Speed_conv_factor={
        'm/s': 1.0,
        'knots': 0.51444444,
    }


    Moment_conv_factor= {
        'N-m':1.0,
        'dyne-cm':1.e-7, 
    }


    
#    cat = raw_input ("Which category would you like to convert (l, p, s, m, t)? ")
#    unit1=raw_input("Which unit would you like to convert from? ")
#    unit2 = raw_input ("Which unit would you like to convert to? ")
#    num1 = int(raw_input ("Enter your value: " ))

    if cat == ("l"):
        conv_dict=Length_conv_factor
    elif cat == ("p"):
        conv_dict = Pressure_Rigidity_conv_factor
    elif cat == ("s"):
        conv_dict = Speed_conv_factor
    elif cat == ("m"):
        conv_dict = Moment_conv_factor

     ## Calculations  

    if cat != ("t"):
        starting_conv_factor = conv_dict.get(unit1)
        standard_num = num1*starting_conv_factor
        final_num = standard_num/conv_dict.get(unit2)
    
    if cat == ("t"):
        if unit1 == 'F':
            if unit2 == 'C':   #  F to C
                final_num = (num1-32)*5/9
            elif unit2 == 'K': #  F to K
                final_num = (num1-32)*5/9+273.15
        elif unit1 == 'C':
            if unit2 == 'F':  #  C to F
                final_num = num1*9/5+32
            elif unit2 == 'K': #  C to K
                final_num = num1+273.15
        elif unit1 == 'K':
            if unit2 == 'F':   #  K to F
                final_num = (num1-273.15)*9/5+32
            elif unit2 == 'C':   #  K to C
                final_num = num1-273.15


    print (final_num)
    return (final_num)


#convertUnits(cat, unit1, unit2, num1)
#convertUnits()   
  

