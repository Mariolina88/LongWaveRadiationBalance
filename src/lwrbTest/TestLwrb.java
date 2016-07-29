package lwrbTest;

import java.net.URISyntaxException;
import java.util.HashMap;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;
import org.jgrasstools.gears.libs.monitor.PrintStreamProgressMonitor;
import org.jgrasstools.hortonmachine.utils.HMTestCase;

import lwrbPointCase.*;



public class TestLwrb extends HMTestCase{

	public TestLwrb() throws Exception {

		String startDate = "2004-06-14 00:00";
		String endDate = "2004-06-16 00:00";
		int timeStepMinutes = 60;
		String fId = "ID";


		PrintStreamProgressMonitor pm = new PrintStreamProgressMonitor(System.out, System.out);


		String inPathToAirT = "resources/input/Taria.csv";
		String inPathToSoilT = "resources/input/Tsuolo.csv";
		String inPathToHumidity = "resources/input/H.csv";
		String inPathToCI = "resources/input/ClearnessIndex.csv";

		String pathToDownwelling= "resources/output/downwelling_model1.csv";
		String pathToUpwelling= "resources/output/upwelling_model1.csv";


		OmsTimeSeriesIteratorReader airTReader = getTimeseriesReader(inPathToAirT, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader soilTReader = getTimeseriesReader(inPathToSoilT, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader humidityReader = getTimeseriesReader(inPathToHumidity, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader CIReader = getTimeseriesReader(inPathToCI, fId, startDate, endDate, timeStepMinutes);
		
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "resources/Input/stations.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;
		
		OmsRasterReader SKYreader = new OmsRasterReader();
		SKYreader.file = "resources/Input/skyview.asc";
		SKYreader.fileNovalue = -9999.0;
		SKYreader.geodataNovalue = Double.NaN;
		SKYreader.process();
		GridCoverage2D skyviewfactor = SKYreader.outRaster;


		OmsTimeSeriesIteratorWriter writer_down = new OmsTimeSeriesIteratorWriter();
		OmsTimeSeriesIteratorWriter writer_up = new OmsTimeSeriesIteratorWriter();



		writer_down.file = pathToDownwelling;
		writer_down.tStart = startDate;
		writer_down.tTimestep = timeStepMinutes;
		writer_down.fileNovalue="-9999";

		writer_up.file = pathToUpwelling;
		writer_up.tStart = startDate;
		writer_up.tTimestep = timeStepMinutes;
		writer_up.fileNovalue="-9999";


		Lwrb lwrb= new Lwrb();


		while( airTReader.doProcess  ) { 

			lwrb.X=-83.89;
			lwrb.Y=1.02;
			lwrb.model="7";
			lwrb.epsilonS=0.98;
			lwrb.A_Cloud=0;
			lwrb.B_Cloud=1;
			
			lwrb.inSkyview = skyviewfactor;
			lwrb.inStations = stationsFC;
			lwrb.fStationsid="cat";

            
			airTReader.nextRecord();	
			HashMap<Integer, double[]> id2ValueMap = airTReader.outData;
			lwrb.inAirTemperatureValues= id2ValueMap;

			soilTReader.nextRecord();
			id2ValueMap = soilTReader.outData;
			lwrb.inSoilTempratureValues = id2ValueMap;

			
			humidityReader.nextRecord();
			id2ValueMap = humidityReader.outData;
			lwrb.inHumidityValues= id2ValueMap;

			CIReader.nextRecord();
			id2ValueMap = CIReader.outData;
			lwrb.inClearnessIndexValues = id2ValueMap;

			lwrb.pm = pm;
			lwrb.process();

			HashMap<Integer, double[]> outHMdown = lwrb.outHMlongwaveDownwelling;
			HashMap<Integer, double[]> outHMup = lwrb.outHMlongwaveUpwelling;



			writer_down.inData = outHMdown;
			writer_down.writeNextLine();

			if (pathToDownwelling != null) {
				writer_down.close();
			}

			writer_up.inData = outHMup;
			writer_up.writeNextLine();

			if (pathToUpwelling != null) {
				writer_up.close();
			}




		}
		airTReader.close();
		soilTReader.close();    
		humidityReader.close();     
		CIReader.close();

	}

	private OmsTimeSeriesIteratorReader getTimeseriesReader( String inPath, String id, String startDate, String endDate,
			int timeStepMinutes ) throws URISyntaxException {
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = inPath;
		reader.idfield = "ID";
		reader.tStart = startDate;
		reader.tTimestep = 60;
		reader.tEnd = endDate;
		reader.fileNovalue = "-9999";
		reader.initProcess();
		return reader;
	}

}