package lwrbTest;

import java.net.URISyntaxException;
import java.util.HashMap;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.io.rasterwriter.OmsRasterWriter;
import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;
import org.jgrasstools.gears.libs.monitor.PrintStreamProgressMonitor;
import org.junit.Test;

import lwrbRasterCase.*;


public class TestLwrbRasterCase{


	GridCoverage2D outDownwellingDataGrid = null;
	GridCoverage2D outUpwellingDataGrid = null;


	@Test
	public void TestLwrb() throws Exception {


		OmsRasterReader skyViewReader = new OmsRasterReader();
		skyViewReader.file = "resources/Input/sky_rid.asc";
		skyViewReader.fileNovalue = -9999.0;
		skyViewReader.geodataNovalue = Double.NaN;
		skyViewReader.process();
		GridCoverage2D skyView = skyViewReader.outRaster;


		OmsRasterReader airTReader = new OmsRasterReader();
		airTReader.file = "resources/Input/airT.asc";
		airTReader.fileNovalue = -9999.0;
		airTReader.geodataNovalue = Double.NaN;
		airTReader.process();
		GridCoverage2D airT = airTReader.outRaster;


		Lwrb lwrb= new Lwrb();

		lwrb.X=5.31;
		lwrb.model="3";
		lwrb.epsilonS=0.98;
		lwrb.A_Cloud=0;
		lwrb.B_Cloud=1;


		lwrb.inSkyviewGrid = skyView;
		lwrb.inAirTempGrid= airT;
		lwrb.inSoilTempGrid= airT;

		lwrb.process();

		outDownwellingDataGrid  = lwrb.outLongwaveDownwellingGrid;
		outUpwellingDataGrid  = lwrb.outLongwaveUpwellingGrid;


		OmsRasterWriter writerDIrectraster = new OmsRasterWriter();
		writerDIrectraster .inRaster = outDownwellingDataGrid;
		writerDIrectraster .file = "resources/Output/downwelling.asc";
		writerDIrectraster.process();

		OmsRasterWriter writerDiffuseraster = new OmsRasterWriter();
		writerDiffuseraster.inRaster = outDownwellingDataGrid;
		writerDiffuseraster.file = "resources/Output/upwelling.asc";
		writerDiffuseraster.process();


	}

}
