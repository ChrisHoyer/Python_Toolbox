

import visa
rm = visa.ResourceManager()
print(rm.list_resources())

#Keysight_PSG = rm.open_resource('GPIB1::20::INSTR')
#zz=Keysight_PSG.query('*IDN?')

#RS_NRP = rm.open_resource('GPIB1::13::INSTR')




#x=RS_NRP.query('READ?')
#RS_NRP.query('FETCh?')
#print(x)

#Betty=rm.open_resource('GPIB1::27::INSTR')
#zz=Betty.query('*IDN?')

#print(zz)
