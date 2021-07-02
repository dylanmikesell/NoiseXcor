## Examples

There are two examples of how to use the NoiseXcor package to process continuous seismic data. 

1. The first is an example of data in SAC format. 
2. The second is an example of how to process data in the MINISEED format.

Within each example are directions on where to download the data files needed to run the examples. 

NOTE: the MINISEED example is currently broken.

### Data folder structure

The raw seismic data need to be organized in such a way that the first script can scan files and populate the database. We chose to follow the different directory structures used in the MSNoise python correlation package.

```
data_structure['BUD']  = "NET/STA/STA.NET.LOC.CHAN.YEAR.DAY"
data_structure['SDS']  = "YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
data_structure['IDDS'] = "YEAR/NET/STA/CHAN.TYPE/DAY/NET.STA.LOC.CHAN.TYPE.YEAR.DAY.HOUR"
data_structure['PDF']  = "YEAR/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
```

We have not yet had data in a format other than BUD, so we have not implemented anything else.

#### Extending to different formats

A user can alter _src/initializeTable.m_ to incorporate new data formats. If someone wants to share other formats, we would be happy to implement them. 
