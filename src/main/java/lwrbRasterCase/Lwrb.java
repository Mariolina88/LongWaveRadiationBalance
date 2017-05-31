/*
 * GNU GPL v3 License
 *
 * Copyright 2015 Marialaura Bancheri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package lwrbRasterCase;


import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.jgrasstools.gears.libs.modules.JGTConstants;

import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;
import java.util.LinkedHashMap;

import javax.media.jai.iterator.RandomIterFactory;
import javax.media.jai.iterator.WritableRandomIter;

import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;
import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.Label;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;
import oms3.annotations.Unit;

import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.utils.RegionMap;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;


import com.vividsolutions.jts.geom.Coordinate;




@Description("The component computes the longwave solar radiation, both upwelling and downwelling.")
@Documentation("")
@Author(name = "Marialaura Bancheri and Giuseppe Formetta", contact = "maryban@hotmail.it")
@Keywords("Hydrology, Radiation, Downwelling , upwelling")
@Label(JGTConstants.HYDROGEOMORPHOLOGY)
@Name("lwrb")
@Status(Status.CERTIFIED)
@License("General Public License Version 3 (GPLv3)")

public class Lwrb extends JGTModel {


	@Description("The map of the interpolated air temperature.")
	@In
	public GridCoverage2D inAirTempGrid;
	
	
	@Description("The map of the soil temprature.")
	@In
	public GridCoverage2D inSoilTempGrid;
	
	@Description("The map of the interpolated humidity.")
	@In
	public GridCoverage2D inHumidityGrid;
	
	
	@Description("The map of the interpolated CI.")
	@In
	public GridCoverage2D inCIGrid;
	
	
	
	
	@Description("Air temperature input value")
	@Unit("°C")
	double airTemperature;


	@Description("Soil temperature input value") 
	@Unit("°C")
	double soilTemperature;


	@Description("Humidity input value") 
	@Unit("%")
	double humidity;

	@Description("Reference humidity")
	private static final double pRH = 0.7;


	@Description("Clearness index input value") 
	@Unit("[0,1]")
	double clearnessIndex;


	@Description("The map of the skyview factor")
	@In
	public GridCoverage2D inSkyviewGrid;
	WritableRaster skyviewfactorWR;

	@Description("X parameter of the literature formulation")
	@In
	public double X; 

	@Description("Y parameter of the literature formulation")
	@In
	public double Y ;

	@Description("Z parameter of the literature formulation")
	@In
	public double Z;


	@Description("Soil emissivity")
	@Unit("-")
	@In
	public double epsilonS;	

	@Description("String containing the number of the model: "
			+ " 1: Angstrom [1918];"
			+ " 2: Brunt's [1932];"
			+ " 3: Swinbank [1963];"
			+ " 4: Idso and Jackson [1969];"
			+ " 5: Brutsaert [1975];"
			+ " 6: Idso [1981];"
			+ " 7: Monteith and Unsworth [1990];"
			+ " 8: Konzelman [1994];"
			+ " 9: Prata [1996];"
			+ " 10: Dilley and O'Brien [1998];"
			+ " 11: To be implemented")
	@In
	public String model;

	@Description("Coefficient to take into account the cloud cover,"
			+ "set equal to 0 for clear sky conditions ")
	@In
	public double A_Cloud;

	@Description("Exponent  to take into account the cloud cover,"
			+ "set equal to 1 for clear sky conditions")
	@In
	public double B_Cloud;

	@Description("It is needed as index of the time step")
	int step;


	@Description("Stefan-Boltzaman costant")
	private static final double ConstBoltz = 5.670373 * Math.pow(10, -8);
	
	
	@Description("The output downwelling map")
	@Out
	public GridCoverage2D outLongwaveDownwellingGrid;
	
	
	@Description("The output upwelling map")
	@Out
	public GridCoverage2D outLongwaveUpwellingGrid;

	Model modelCS;
	WritableRaster airTemperatureMap;
	WritableRaster soilTemperatureMap;
	WritableRaster humidityMap;
	WritableRaster CIMap;
	
	int cols;
	int rows;
	double dx;
	RegionMap regionMap;
	
	@Description("the linked HashMap with the coordinate of the stations")
	LinkedHashMap<Integer, Coordinate> stationCoordinates;



	/**
	 * Process.
	 *
	 * @throws Exception the exception
	 */
	@Execute
	public void process() throws Exception {

		if(step==0){
			
			airTemperatureMap=mapsTransform(inAirTempGrid);
			soilTemperatureMap=mapsTransform(inSoilTempGrid);
			if (CIMap!= null) CIMap=mapsTransform(inCIGrid);
			if (humidityMap!= null) humidityMap=mapsTransform(inHumidityGrid);
			skyviewfactorWR =mapsTransform(inSkyviewGrid);

			// get the dimension of the maps
			regionMap = CoverageUtilities.getRegionParamsFromGridCoverage(inSkyviewGrid);
			cols = regionMap.getCols();
			rows = regionMap.getRows();
			dx=regionMap.getXres();
		}
		
		WritableRaster outDownwellingWritableRaster = CoverageUtilities.createDoubleWritableRaster(cols, rows, null, null, null);
		WritableRaster outUpwellingWritableRaster =CoverageUtilities.createDoubleWritableRaster(cols, rows, null, null, null); 
	

		WritableRandomIter downwellingIter = RandomIterFactory.createWritable(outDownwellingWritableRaster, null);       
		WritableRandomIter upwellingIter = RandomIterFactory.createWritable(outUpwellingWritableRaster, null);		



		// get the geometry of the maps and the coordinates of the stations
		GridGeometry2D inSkyGridGeo = inSkyviewGrid.getGridGeometry();
		stationCoordinates = getCoordinate(inSkyGridGeo);
		
		
		// iterate over the entire domain and compute for each pixel the SWE
		for( int r = 1; r < rows - 1; r++ ) {
			for( int c = 1; c < cols - 1; c++ ) {

			airTemperature=airTemperatureMap.getSampleDouble(c, r, 0);

			soilTemperature = soilTemperatureMap.getSampleDouble(c, r, 0);

			humidity= pRH;
			if (humidityMap!= null)
				humidity = humidityMap.getSampleDouble(c, r, 0);
			if(isNovalue(humidity)) humidity=pRH;

			clearnessIndex = 1;
			if (CIMap != null) clearnessIndex = CIMap.getSampleDouble(c, r, 0);
			if (isNovalue(clearnessIndex )) clearnessIndex = 1;

			double skyviewvalue=skyviewfactorWR.getSampleDouble(c, r, 0);

			/**Computation of the downwelling, upwelling and longwave:
			 * if there is no value in the input data, there will be no value also in
			 * the output*/
			
			double upwelling=(isNovalue(soilTemperature))? Double.NaN:computeUpwelling(soilTemperature);
			upwelling=(upwelling<0)? Double.NaN:upwelling;
			upwelling=(upwelling>2000)? Double.NaN:upwelling;
			
			
			double downwellingALLSKY=(isNovalue(airTemperature))? Double.NaN:
				computeDownwelling(model,airTemperature,humidity/100,skyviewvalue,upwelling);
			
			downwellingALLSKY=(downwellingALLSKY<0)? Double.NaN:downwellingALLSKY;
			downwellingALLSKY=(downwellingALLSKY>2000)? Double.NaN:downwellingALLSKY;
			
			downwellingIter.setSample(c, r, 0,downwellingALLSKY);
			upwellingIter.setSample(c, r, 0,upwelling);
			}
		}
		
		CoverageUtilities.setNovalueBorder(outDownwellingWritableRaster);
		CoverageUtilities.setNovalueBorder(outUpwellingWritableRaster);	

		

		outLongwaveDownwellingGrid = CoverageUtilities.buildCoverage("Downwelling", outDownwellingWritableRaster, 
				regionMap, inSkyviewGrid.getCoordinateReferenceSystem());
		outLongwaveUpwellingGrid= CoverageUtilities.buildCoverage("Upwelling", outUpwellingWritableRaster, 
				regionMap, inSkyviewGrid.getCoordinateReferenceSystem());

		step++;
	}

	
	/**
	 * Gets the coordinate of each pixel of the given map.
	 *
	 * @param GridGeometry2D grid is the map 
	 * @return the coordinate of each point
	 */
	private LinkedHashMap<Integer, Coordinate> getCoordinate(GridGeometry2D grid) {
		LinkedHashMap<Integer, Coordinate> out = new LinkedHashMap<Integer, Coordinate>();
		int count = 0;
		RegionMap regionMap = CoverageUtilities.gridGeometry2RegionParamsMap(grid);
		double cols = regionMap.getCols();
		double rows = regionMap.getRows();
		double south = regionMap.getSouth();
		double west = regionMap.getWest();
		double xres = regionMap.getXres();
		double yres = regionMap.getYres();
		double northing = south;
		double easting = west;
		for (int i = 0; i < cols; i++) {
			easting = easting + xres;
			for (int j = 0; j < rows; j++) {
				northing = northing + yres;
				Coordinate coordinate = new Coordinate();
				coordinate.x = west + i * xres;
				coordinate.y = south + j * yres;
				out.put(count, coordinate);
				count++;
			}
		}

		return out;
	}


	/**
	 * Compute upwelling longwave radiation .
	 *
	 * @param soilTemperature: the soil temperature input
	 * @return the double value of the upwelling
	 */
	private double computeUpwelling( double soilTemperature){

		/**compute the upwelling*/
		return epsilonS * ConstBoltz * Math.pow(soilTemperature+ 273.15, 4);
	}

	/**
	 * Compute downwelling longwave radiation.
	 *
	 * @param model: the string containing the number of the model
	 * @param airTemperature:  the air temperature input
	 * @param humidity: the humidity input
	 * @param clearnessIndex: the clearness index input
	 * @return the double value of the all sky downwelling
	 */
	private double computeDownwelling(String model,double airTemperature, 
			double humidity, double skyviewvalue, double upwelling){

		/**e is the screen-level water-vapor pressure in kPa*/
		double e = humidity *6.11 * Math.pow(10, (7.5 * airTemperature) / (237.3 + airTemperature)) / 10;

		/**compute the clear sky emissivity*/
		modelCS=SimpleModelFactory.createModel(model,X,Y,Z,airTemperature+ 273.15,e);
		double epsilonCS=modelCS.epsilonCSValues();

		/**compute the downwelling in clear sky conditions*/
		double downwellingCS=epsilonCS* ConstBoltz* Math.pow(airTemperature+ 273.15, 4);
		
		/**correct downwelling clear sky for sloping terrain*/
		downwellingCS=downwellingCS*skyviewvalue+upwelling*(1-skyviewvalue);

		/**compute the cloudness index*/
		double cloudnessIndex = 1 + A_Cloud* Math.pow((1-clearnessIndex), B_Cloud);

		/**compute the downwelling in all-sky conditions*/
		return downwellingCS * cloudnessIndex;

	}

	
	/**
	 * Maps reader transform the GrifCoverage2D in to the writable raster and
	 * replace the -9999.0 value with no value.
	 *
	 * @param inValues: the input map values
	 * @return the writable raster of the given map
	 */
	private WritableRaster mapsTransform  ( GridCoverage2D inValues){	
		RenderedImage inValuesRenderedImage = inValues.getRenderedImage();
		WritableRaster inValuesWR = CoverageUtilities.replaceNovalue(inValuesRenderedImage, -9999.0);
		inValuesRenderedImage = null;
		return inValuesWR;
	}


}