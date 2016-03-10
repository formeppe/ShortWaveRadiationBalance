/*
 * This file is part of JGrasstools (http://www.jgrasstools.org)
 * (C) HydroloGIS - www.hydrologis.com 
 * 
 * JGrasstools is free software: you can redistribute it and/or modify
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
package swrbRasterCase;

import static org.jgrasstools.gears.libs.modules.ModelsEngine.calcInverseSunVector;
import static org.jgrasstools.gears.libs.modules.ModelsEngine.calcNormalSunVector;
import static org.jgrasstools.gears.libs.modules.ModelsEngine.calculateFactor;
import static org.jgrasstools.gears.libs.modules.ModelsEngine.scalarProduct;

import org.jgrasstools.gears.libs.modules.JGTConstants;

import java.awt.image.RenderedImage;
import java.awt.image.SampleModel;
import java.awt.image.WritableRaster;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

import static org.jgrasstools.gears.libs.modules.JGTConstants.doubleNovalue;
import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;

import javax.media.jai.RasterFactory;
import javax.media.jai.iterator.RandomIter;
import javax.media.jai.iterator.RandomIterFactory;
import javax.media.jai.iterator.WritableRandomIter;

import oms3.annotations.Author;
import oms3.annotations.Bibliography;
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

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.geometry.DirectPosition2D;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.jgrasstools.gears.io.rasterwriter.OmsRasterWriter;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.libs.monitor.LogProgressMonitor;
import org.jgrasstools.gears.utils.CrsUtilities;
import org.jgrasstools.gears.utils.RegionMap;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.gears.utils.geometry.GeometryUtilities;
import org.joda.time.DateTime;
import org.joda.time.DateTimeZone;
import org.joda.time.format.DateTimeFormat;
import org.joda.time.format.DateTimeFormatter;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.geometry.DirectPosition;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Point;

@Description("Calculate the amount of beam and diffuse shortwave radiation .")
@Documentation("")
@Author(name = "Marialaura Bancheri, Giuseppe Formetta, Daniele Andreis and Riccardo Rigon", contact = "maryban@hotmail.it")
@Keywords("Hydrology, Radiation, SkyviewFactor, Hillshade")
@Bibliography("Formetta (2013)")
@Label(JGTConstants.HYDROGEOMORPHOLOGY)
@Name("shortradbal")
@Status(Status.CERTIFIED)
@License("General Public License Version 3 (GPLv3)")

public class ShortwaveRadiationBalanceRasterCase extends JGTModel {

	@Description("The map of the interpolated temperature.")
	@In
	public GridCoverage2D inTempGrid;

	@Description("The double value of the  temperature")
	double temperature;

	@Description("The map of the the interpolated humidity.")
	@In
	public GridCoverage2D inHumidityGrid;

	@Description("The double value of the humidity")
	double humidity;

	@Description("The digital elevation model.")
	@In
	public GridCoverage2D inDem;
	WritableRaster demWR;
	
	@Description("The map of the skyview factor")
	@In
	public GridCoverage2D inSkyview;
	WritableRaster skyviewfactorWR;
	
	@Description("The timeStep allows to chose between the hourly time step"
			+ " or the daily time step. It could be: "
			+ " Hourly or Daily")
	@In
	public String timeStep;

	@Description("It is needed to iterate on the date")
	int step;

	@Description("The first day of the simulation.")
	@In
	public String tStartDate;

	DateTimeFormatter formatter = DateTimeFormat.forPattern("yyyy-MM-dd HH:mm").withZone(DateTimeZone.UTC);

	@Description("Final target CRS")
	CoordinateReferenceSystem targetCRS = DefaultGeographicCRS.WGS84;

	@Description("Ozone layer thickness in cm")
	@In
	public double pCmO3;
	// pCmO3 = 0.4-0.6;

	@Description("Default relative humidity value")
	public double pRH = 0.7;

	@Description(" For aerosol attenuation (5 < vis < 180 Km) [km].")
	@In
	@Unit ("km")
	public double pVisibility;
	//pVisibility = 60-80;

	@Description("The soil albedo.")
	@In
	public double pAlphag;
	//pAlphag = 0.9;

	@Description("The solar constant")
	private static final double SOLARCTE = 1370.0;
	// double SOLARCTE = 1360.0;

	@Description("The atmospheric pressure")
	private static final double ATM = 1013.25;

	@Description("The declination of the sun, diven the day")
	double delta;

	@Description("the linked HashMap with the coordinate of the stations")
	LinkedHashMap<Integer, Coordinate> stationCoordinates;

	@Description("List of the latitudes of the station ")
	ArrayList <Double> latitudeStation= new ArrayList <Double>();
	
	@Description("relative air mass")
	double ma;

	@Description("trasmittance function for the Reyleigh scattering")
	double tau_r;

	@Description("trasmittance by ozone")
	double tau_o;

	@Description("trasmiitance by gases")
	double tau_g;

	@Description("trasmittance by aerosol")
	double tau_a;

	@Description("trasmittance by water vapor")
	double tau_w;

	@Description("direct normal irradiance")
	double In;

	@Description("the hour of the consdiered day")
	double hour;

	@Description("the sunrise in the considered day")
	double sunrise;

	@Description("the sunrise in the considered day")
	double sunset;
	
	@Description("The output direct radiation map")
	@Out
	public GridCoverage2D outDirectGrid;

	@Description("The output diffuse radiation map")
	@Out
	public GridCoverage2D outDiffuseGrid;
	
	@Description("The output top atmosphere map")
	@Out
	public GridCoverage2D outTopATMGrid;
	
	

	@Execute
	public void process() throws Exception {

		// This 2 operations allow to define if we are working with daily or hourly time step
		// if we are working with Daily time step, every time it adds to the start date a day
		// otherwise it adds an hour, "step increments at the end of the process
		// the actual date is needed to compute the actual energy index	
		DateTime startDateTime = formatter.parseDateTime(tStartDate);
		DateTime date=(timeStep.equals("Daily"))?startDateTime.plusDays(step):startDateTime.plusHours(step);

		// computing the reference system of the input DEM
		CoordinateReferenceSystem sourceCRS = inDem.getCoordinateReferenceSystem2D();

		// transform the GrifCoverage2D maps into writable rasters
		WritableRaster temperatureMap=mapsReader(inTempGrid);	
		WritableRaster humidityMap=mapsReader(inHumidityGrid);
		demWR=mapsReader(inDem);
		
		// get the dimension of the maps
		RegionMap regionMap = CoverageUtilities.getRegionParamsFromGridCoverage(inDem);
		int cols = regionMap.getCols();
		int rows = regionMap.getRows();
		double dx=regionMap.getXres();
		

		WritableRaster outDiffuseWritableRaster = CoverageUtilities.createDoubleWritableRaster(cols, rows, null, null, null);
		WritableRaster outDirectWritableRaster =CoverageUtilities.createDoubleWritableRaster(cols, rows, null, null, null); 
		WritableRaster outTopATMWritableRaster =CoverageUtilities.createDoubleWritableRaster(cols, rows, null, null, null); 

		WritableRandomIter diffuseIter = RandomIterFactory.createWritable(outDiffuseWritableRaster, null);       
		WritableRandomIter directIter = RandomIterFactory.createWritable(outDirectWritableRaster, null);		
		WritableRandomIter topIter = RandomIterFactory.createWritable(outTopATMWritableRaster, null);       


		// get the geometry of the maps and the coordinates of the stations
		GridGeometry2D inDemGridGeo = inDem.getGridGeometry();
		stationCoordinates = getCoordinate(inDemGridGeo);

		// iterate over the entire domain and compute for each pixel the SWE
		for( int r = 1; r < rows - 1; r++ ) {
			for( int c = 1; c < cols - 1; c++ ) {
				int k=0;

				// get the exact value of the variable in the pixel i, j 
				temperature=temperatureMap.getSampleDouble(c, r, 0);
				humidity=humidityMap.getSampleDouble(c, r, 0);
				if(isNovalue(humidity)) humidity=pRH;

				// get the coordinate of the given pixel
				Coordinate coordinate = (Coordinate) stationCoordinates.get(k);

				//compute the latitude, after a reprojection in WGS 84
				Point [] idPoint=getPoint(coordinate,sourceCRS, targetCRS);
				latitudeStation.add(Math.toRadians(idPoint[0].getY()));

				// calculating the sun vector
				double sunVector[] = calcSunVector(latitudeStation.get(k), getHourAngle(date,latitudeStation.get(k)));

				// calculate the inverse of the sun vector
				double[] inverseSunVector = calcInverseSunVector(sunVector);

				// calculate the normal of the sun vector according to Corripio (2003)
				double[] normalSunVector = calcNormalSunVector(sunVector);

				//E0 is the correction factor related to Earth’s orbit eccentricity computed according to Spencer (1971):
				double E0=computeE0(date);
				
				//evaluate the shadow map
				WritableRaster shadowWR = calculateFactor(rows,cols,  sunVector, inverseSunVector, normalSunVector, demWR, dx);

				// compute the vector normal to a grid cell surface.
				WritableRaster normalWR = normalVector(demWR, dx);
				

				// calculate the direct radiation, during the daylight
				double direct= (hour > (sunrise) && hour < (sunset))?
						calcDirectRadiation(c, r, demWR, shadowWR, sunVector, normalWR,E0,temperature, humidity):0;

				//calculate the diffuse radiation, during the daylight
				double diffuse =(hour > (sunrise) && hour < (sunset))?
								calcDiffuseRadiation(sunVector, E0,c, r):0;

				// calculate the radiation at the top of the atmosphere, during the daylight
				double topATM=(hour > (sunrise) && hour < (sunset))?
								calcTopAtmosphere(E0, sunVector[2]):0;
								
				double convert=(timeStep.equals("Daily"))?277.7:11.57;
								
				diffuseIter.setSample(c, r, 0,direct/convert);
				directIter.setSample(c, r, 0, diffuse/convert);
				topIter.setSample(c, r, 0, topATM/convert);

				// the index k is for the loop over the list
				k++;


			}
		}

		CoverageUtilities.setNovalueBorder(outDirectWritableRaster);
		CoverageUtilities.setNovalueBorder(outDiffuseWritableRaster);	
		CoverageUtilities.setNovalueBorder(outTopATMWritableRaster);

		outDirectGrid = CoverageUtilities.buildCoverage("Direct", outDirectWritableRaster, 
				regionMap, inDem.getCoordinateReferenceSystem());
		outDiffuseGrid= CoverageUtilities.buildCoverage("Diffuse", outDiffuseWritableRaster, 
				regionMap, inDem.getCoordinateReferenceSystem());
		
		outTopATMGrid= CoverageUtilities.buildCoverage("topATM", outTopATMWritableRaster, 
				regionMap, inDem.getCoordinateReferenceSystem());

		// upgrade the step for the new date
		step++;	

	}

	/**
	 * Maps reader transform the GrifCoverage2D in to the writable raster and
	 * replace the -9999.0 value with no value.
	 *
	 * @param inValues: the input map values
	 * @return the writable raster of the given map
	 */
	private WritableRaster mapsReader ( GridCoverage2D inValues){	
		RenderedImage inValuesRenderedImage = inValues.getRenderedImage();
		WritableRaster inValuesWR = CoverageUtilities.replaceNovalue(inValuesRenderedImage, -9999.0);
		inValuesRenderedImage = null;
		return inValuesWR;
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
	 * Gets the point.
	 *
	 * @param coordinate the coordinate
	 * @param sourceCRS is the source crs
	 * @param targetCRS the target crs
	 * @return the point
	 * @throws Exception the exception
	 */
	private Point[] getPoint(Coordinate coordinate, CoordinateReferenceSystem sourceCRS, CoordinateReferenceSystem targetCRS) 
			throws Exception{
		Point[] point = new Point[] { GeometryUtilities.gf().createPoint(coordinate) };
		CrsUtilities.reproject(sourceCRS, targetCRS, point);
		return point;
	}

	/**
	 * Compute the correction factor related to Earth’s orbit eccentricity.
	 *
	 * @param date is the current date
	 * @return the double value of E0
	 */


	private double computeE0(DateTime date) {
		// k is the day angle in radiant 
		double k = 2 * Math.PI * (date.getDayOfMonth() - 1.0) / 365.0;
		return 1.00011 + 0.034221 * Math.cos(k) + 0.00128
				* Math.sin(k) + 0.000719 * Math.cos(2 * k) + 0.000077
				* Math.sin(2 * k);
	}
	
	/**
	 * getHourAngle is the value of the hour angle at a given time and latitude (Corripio (2003))
	 * 
	 *
	 * @param date is the current date
	 * @param latitude is the latitude of the station
	 * @return the double value of the hour angle
	 */
	private double getHourAngle(DateTime date, double latitude) {
		int day = date.getDayOfYear();	
		hour=(double)date.getMillisOfDay() / (1000 * 60 * 60);

		// (360 / 365.25) * (day - 79.436) is the number of the day 
		double dayangb = Math.toRadians((360 / 365.25) * (day - 79.436));

		// Evaluate the declination of the sun.
		delta = Math.toRadians(.3723 + 23.2567 * Math.sin(dayangb) - .758
				* Math.cos(dayangb) + .1149 * Math.sin(2 * dayangb) + .3656
				* Math.cos(2 * dayangb) - .1712 * Math.sin(3 * dayangb) + .0201
				* Math.cos(3 * dayangb));

		// ss is the absolute value of the hour angle at sunrise or sunset
		double ss = Math.acos(-Math.tan(delta) * Math.tan(latitude));
		sunrise = 12 * (1.0 - ss / Math.PI);
		sunset = 12 * (1.0 + ss / Math.PI);

		if (hour > (sunrise) && hour < (sunset) & (hour - sunrise) < 0.01) hour = hour + 0.1;
		if (hour > (sunrise) && hour < (sunset) & (sunset - hour) < 0.01) hour = hour - 0.1;

		//the hour angle is zero at noon and has the following value in radians at any time t
		// given in hours and decimal fraction:
		double hourangle=(hour/ 12.0 - 1.0) * Math.PI;
		return hourangle;

	}
	
	/**
	 * calcSunVector compute the vector vector in the direction of the Sun (Corripio (2003))
	 *
	 * @param latitude is the latitude of the station 
	 * @param the hour angle  
	 * @return the sun vector 
	 */
	protected double[] calcSunVector(double latitude, double hourAngle) {
		double sunVector[] = new double[3];
		sunVector[0] = -Math.sin(hourAngle) * Math.cos(delta);
		sunVector[1] = Math.sin(latitude) * Math.cos(hourAngle) * Math.cos(delta)
				- Math.cos(latitude) * Math.sin(delta);
		sunVector[2] = Math.cos(latitude) * Math.cos(hourAngle) * Math.cos(delta)
				+ Math.sin(latitude) * Math.sin(delta);
		return sunVector;

	}
	
	/**
	 * normalVector compute the vector normal to a grid cell surface, according to Corripio (2003)
	 *
	 * @param demWR is the Writable raster of DEM 
	 * @param res is the resolution of the DEM
	 * @return the normal vector for each cell
	 */
	protected WritableRaster normalVector(WritableRaster demWR, double res) {

		int minX = demWR.getMinX();
		int minY = demWR.getMinY();
		int rows = demWR.getHeight();
		int cols = demWR.getWidth();

		RandomIter pitIter = RandomIterFactory.create(demWR, null);
		/*
		 * Initialize the image of the normal vector in the central point of the
		 * cells, which have 3 components (X;Y;Z), so the Image have 3 bands..
		 */
		SampleModel sm = RasterFactory.createBandedSampleModel(5, cols, rows, 3);
		WritableRaster tmpNormalVectorWR = CoverageUtilities .createDoubleWritableRaster(cols, rows, null, sm, 0.0);
		WritableRandomIter tmpNormalIter = RandomIterFactory.createWritable( tmpNormalVectorWR, null);
		/*
		 * apply the corripio's formula 
		 */
		for (int j = minY; j < minX + rows - 1; j++) {
			for (int i = minX; i < minX + cols - 1; i++) {
				double zij = pitIter.getSampleDouble(i, j, 0);
				double zidxj = pitIter.getSampleDouble(i + 1, j, 0);
				double zijdy = pitIter.getSampleDouble(i, j + 1, 0);
				double zidxjdy = pitIter.getSampleDouble(i + 1, j + 1, 0);
				double firstComponent = res * (zij - zidxj + zijdy - zidxjdy);
				double secondComponent = res * (zij + zidxj - zijdy - zidxjdy);
				double thirthComponent = 2 * (res * res);
				double den = Math.sqrt(firstComponent * firstComponent
						+ secondComponent * secondComponent + thirthComponent
						* thirthComponent);
				tmpNormalIter.setPixel(i, j, new double[] {
						firstComponent / den, secondComponent / den,
						thirthComponent / den });

			}
		}
		pitIter.done();

		return tmpNormalVectorWR;

	}

	/**
	 * calcDirectRadiation calculates the direct radiation according to Corripio (2002)
	 *
	 * @param i is the the index of the column of the station in the DEM
	 * @param j is the the index of the row of the station in the DEM
	 * @param demWR is the DEM
	 * @param shadowWR is the map with the shadow 
	 * @param sunVector is the sun vector
	 * @param normalWR is the the raster with the vector normal to a pixel
	 * @param E0 is the correction of the eccentricity
	 * @param temperture is the air temperature
	 * @param himidity is the relative humidity
	 * @return the double value of the direct radiation
	 */

	private double calcDirectRadiation(int i, int j, WritableRaster demWR, WritableRaster shadowWR, double[] sunVector,
			WritableRaster normalWR,double E0, double temperature, double humidity) throws IOException {

		// zenith angle
		double zenith = Math.acos(sunVector[2]);

		//mr [–] relative optical air mass:
		double mr = 1.0 / (sunVector[2] + 0.15 * Math.pow( (93.885 - (zenith * (180 / (2*Math.PI)))), (-1.253)));

		// altitude of the station
		double z = demWR.getSampleDouble(i, j, 0);

		// local atmospheric pressure
		double pressure = ATM * Math.exp(-0.0001184 * z);

		// relative air mass
		ma = mr * pressure / ATM;

		// The transmittance functions for Rayleigh scattering
		tau_r = Math.exp((-.09030 * Math.pow(ma, 0.84))
				* (1.0 + ma - Math.pow(ma, 1.01)));

		//transform the temperature in Kelvin
		temperature  = temperature + 273.0;

		// evaluate the saturated valor pressure
		double saturatedVaporPressure = Math.exp(26.23 - 5416.0 / temperature );

		// the precipitable water in cm calculated according to Prata (1996)
		double w = 0.493 * (humidity / 100) * saturatedVaporPressure / temperature ;

		// Transmittance by water vapour
		tau_w = 1.0 - 2.4959 * w * mr / (Math.pow(1.0 + 79.034 * w * mr, 0.6828) + 6.385 * w * mr);

		// The transmittance by ozone 
		tau_o = 1.0 - ((0.1611 * pCmO3 * mr * Math.pow(1.0 + 139.48 * pCmO3 * mr,-0.3035)) 
				- (0.002715 * pCmO3 * mr / (1.0 + 0.044 * pCmO3 * mr + 0.0003 * Math.pow(pCmO3 * mr, 2))));

		// Transmittance by uniformly mixed gases
		tau_g = Math.exp(-0.0127 * Math.pow(ma, 0.26));


		// The transmittance by aerosols
		tau_a = Math.pow((0.97 - 1.265 * Math.pow(pVisibility,(-0.66))), Math.pow(ma, 0.9));

		// correction factor [m] for increased trasmittance with elevation z[m] according to Corripio (2002)
		double beta_s = (z <= 3000)?2.2 * Math.pow(10, -5) * z:2.2 * Math.pow(10, -5) * 3000;


		// cosin: incidence angle according to Corripio (2003)
		double cos_inc = scalarProduct(sunVector, normalWR.getPixel(i, j, new double[3]));
		if (cos_inc < 0) cos_inc = 0;

		// Direct radiation under cloudless sky incident on arbitrary tilted
		// surfaces (by inclusion of cosinc)
		In=0.9571*SOLARCTE*E0*(tau_r * tau_o * tau_g * tau_w * tau_a + beta_s);

		double S_incident=In* cos_inc * shadowWR.getSampleDouble(i, j, 0);

		return S_incident=(S_incident<=0)?doubleNovalue:S_incident;
	}


	/**
	 * calcDiffuseRadiation calculates the diffuse radiation according to Corripio (2002)
	 *
	 * @param i is the the index of the column of the station in the DEM
	 * @param j is the the index of the row of the station in the DEM 
	 * @param sunVector is the sun vector
	 * @param E0 is the correction of the eccentricity
	 * @return the double value of the diffuse radiation
	 */
	private double calcDiffuseRadiation(double [] sunVector,double E0, int i, int j){

		//single-scattering albedo fraction of incident energy scattered to total attenuation by aerosol
		double omega0 = 0.9;

		//trasmittance of direct radiation due to aerosol absorbance 
		double tau_aa = 1.0 - (1.0 - omega0)* (1.0 - ma + Math.pow(ma, 1.06)) * (1 - tau_a);
		// ////////////////////////////////////////////////////

		//  Rayleigh scattered diffunce irradiance
		double I_dr = 0.79 * E0*SOLARCTE* sunVector[2] * (1.0 - tau_r)* (tau_o * tau_g * tau_w * tau_aa) * 0.5/ (1.0 - ma + Math.pow(ma, 1.02));

		// The aerosol-scattered diffuse irradiance 
		double FC = 0.74;
		double I_da = 0.79 * E0*SOLARCTE* sunVector[2] * (tau_o * tau_g * tau_w * tau_aa) * FC
				* (1.0 - (tau_a / tau_aa)) / ((1 - ma + Math.pow(ma, 1.02)));

		// The atmospheric albedo is computed as
		double alpha_a = 0.0685 + (1.0 - FC) * (1.0 - (tau_a / tau_aa));

		// the diffuse irradiance from multiple reflection between the earth and the atmosphere
		double I_dm = (In*sunVector[2]+ I_dr + I_da) * alpha_a * pAlphag / (1.0 - pAlphag * alpha_a);

		double diffuse = (I_dr + I_da + I_dm)* skyviewfactorWR.getSampleDouble(i, j, 0);

		return diffuse=(diffuse<0)?doubleNovalue:diffuse;

	}

	/**
	 * calcTopAtmosphere calculates the radiation at the top of the atmosphere according to Corripio (2002)
	 *
	 * @param sunVector is the sun vector
	 * @param E0 is the correction of the eccentricity
	 * @return the double value of the top atmosphere radiation
	 */
	private double calcTopAtmosphere(double E0, double sunVector){
		double topATM=E0 * SOLARCTE * Math.cos(Math.acos(sunVector));
		return topATM=(topATM<0)?0:topATM;	
	}

}