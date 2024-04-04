# ukesm_to_jules

My code is here, it calls Karina’s script:
/home/h06/eroberts/workspace/UM_TO_JULES/um_to_jules_karina.py

It uses a local copy of her script which I’ve modified a bit.

https://code.metoffice.gov.uk/trac/UKESM/wiki/AtmosDiagsForcingOceanJules



It was Mule that I’d used previously to edit a year or some other aspect of header info (on xcslr1 in /home/d04/hadsl/UKESM1 ) . Apparently you can manipulate fields with it, the documentation is here, with example scripts:

https://code.metoffice.gov.uk/doc/um/mule/latest/index.html
https://code.metoffice.gov.uk/doc/um/mule/latest/user_guide.html
https://code.metoffice.gov.uk/doc/um/mule/latest/examples.html



import mule

ff1 = mule.FieldsFile.from_file("InputDump")  #full UM dump to be modified
ff2 = mule.FieldsFile.from_file("InputAncil") #single field ancil, we’ve already made from a JULES dump

ff_out = ff1.copy()

for field_1 in ff1.fields:     #loop over all fields in dump
    if field_1.lbuser4 == 466: #if stash code is for first soil carbon pool
       ff_out.fields.remove(field_1) #remove the old soil carbon data
       break

ff_out.fields.append(ff2) #add the new soil carbon data

ff_out.to_file("OutputFile") #save the modified dump
