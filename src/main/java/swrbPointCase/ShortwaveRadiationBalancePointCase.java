/*
 * GNU GPL v3 License
 *
 * Copyright 2022 Giuseppe Formetta
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
package swrbPointCase;

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
import java.util.Calendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;

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
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.feature.SchemaException;
import org.geotools.geometry.DirectPosition2D;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.utils.CrsUtilities;
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
@Author(name = "Giuseppe Formetta, Marialaura Bancheri, Daniele Andreis and Riccardo Rigon", contact = "giuseppe.formetta@unitn.it")
@Keywords("Hydrology, Radiation, SkyviewFactor, Hillshade")
@Bibliography("Formetta (2013)")
@Label(JGTConstants.HYDROGEOMORPHOLOGY)
@Name("shortradbal")
@Status(Status.CERTIFIED)
@License("General Public License Version 3 (GPLv3)")
public class ShortwaveRadiationBalancePointCase extends JGTModel {

	@Description("The Hashmap with the time series of the temperature values")
	@In
	@Unit("°C")
	public HashMap<Integer, double[]> inTemperatureValues;

	@Description("The double value of the  temperature, once read from the HashMap")
	double temperature;

	@Description("The Hashmap with the time series of the humidity values")
	@In
	@Unit("%")
	public HashMap<Integer, double[]> inHumidityValues;

	@Description("The double value of the  humidity, once read from the HashMap")
	double humidity;

	@Description("The map of the Digital Elevation Model")
	@In
	public GridCoverage2D inDem;
	WritableRaster demWR;

	@Description("The map of the skyview factor")
	@In
	public GridCoverage2D inSkyview;
	WritableRaster skyviewfactorWR;

	@Description("The shape file with the station measuremnts")
	@In
	public SimpleFeatureCollection inStations;

	@Description("The name of the field containing the ID of the station in the shape file")
	@In
	public String fStationsid;

	@Description(" The vetor containing the id of the station")
	Object[] idStations;

	@Description("the linked HashMap with the coordinate of the stations")
	LinkedHashMap<Integer, Coordinate> stationCoordinates;

	@Description("doHourly allows to chose between the hourly time step" + " or the daily time step. It could be: "
			+ " Hourly--> true or Daily-->false")
	@In
	public boolean doHourly;

	@Description("It is needed to iterate on the date")
	int step;

	@Description("The first day of the simulation.")
	@In
	public String tStartDate;

	@Description("Ozone layer thickness in cm")
	@In
	public double pCmO3;
	// pCmO3 = 0.4-0.6;

	@Description("Default relative humidity value")
	public double pRH = 0.7;

	@Description(" For aerosol attenuation (5 < vis < 180 Km) [km].")
	@In
	@Unit("km")
	public double pVisibility;
	// pVisibility = 60-80;

	@Description("The soil albedo.")
	@In
	public double pAlphag;
	// pAlphag = 0.9;

	@Description("The temperature default value in case of missing data.")
	@In
	@Unit("C")
	public double defaultTemp = 15.0;

	@Description("The solar constant")
	private static final double SOLARCTE = 1370.0;
	// double SOLARCTE = 1360.0;

	@Description("The atmospheric pressure")
	private static final double ATM = 1013.25;

	@Description("The declination of the sun, diven the day")
	double delta;

	@Description("List of the indeces of the columns of the station in the map")
	ArrayList<Integer> columnStation = new ArrayList<Integer>();

	@Description("List of the indeces of the rows of the station in the map")
	ArrayList<Integer> rowStation = new ArrayList<Integer>();

	@Description("List of the latitudes of the station ")
	ArrayList<Double> latitudeStation = new ArrayList<Double>();

	@Description("List of the latitudes of the station ")
	ArrayList<Double> longitudeStation = new ArrayList<Double>();

	@Description("Final target CRS")
	CoordinateReferenceSystem targetCRS = DefaultGeographicCRS.WGS84;

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

	@Description("the hour of the consdired day")
	double hour;

	@Description("the sunrise in the considered day")
	double sunrise;

	@Description("the sunrise in the considered day")
	double sunset;

	@Description("the output hashmap withe the direct radiation")
	@Out
	public HashMap<Integer, double[]> outHMdirect = new HashMap<Integer, double[]>();

	@Description("the output hashmap withe the diffuse radiation")
	@Out
	public HashMap<Integer, double[]> outHMdiffuse = new HashMap<Integer, double[]>();

	@Description("the output hashmap withe the top atmosphere radiation")
	@Out
	public HashMap<Integer, double[]> outHMtopatm = new HashMap<Integer, double[]>();

	@Description("the output hashmap withe the top atmosphere radiation")
	@Out
	public HashMap<Integer, double[]> outHMtotal = new HashMap<Integer, double[]>();

	DateTimeFormatter formatter = DateTimeFormat.forPattern("yyyy-MM-dd HH:mm").withZone(DateTimeZone.UTC);

	WritableRaster normalWR;
	int height;
	int width;
	double dx;

	@Execute
	public void process() throws Exception {

		// This 2 operations allow to define if we are working with daily or hourly time
		// step
		// if we are working with Daily time step, every time it adds to the start date
		// a day
		// otherwise it adds an hour, "step" increments at the end of the process
		// the actual date is needed to compute the actual sunrise and sunset
		DateTime startDateTime = formatter.parseDateTime(tStartDate);
		DateTime date = (doHourly == false) ? startDateTime.plusDays(step)
				: startDateTime.plusHours(step).plusMinutes(30);

		// from pixel coordinates (in coverage image) to geographic coordinates (in
		// coverage CRS)
		MathTransform transf = inDem.getGridGeometry().getCRSToGrid2D();

		// computing the reference system of the input DEM
		CoordinateReferenceSystem sourceCRS = inDem.getCoordinateReferenceSystem2D();

		if (step == 0) {
			// transform the GrifCoverage2D maps into writable rasters
			demWR = mapsTransform(inDem);
			skyviewfactorWR = mapsTransform(inSkyview);

			// starting from the shp file containing the stations, get the coordinate
			// of each station
			stationCoordinates = getCoordinate(inStations, fStationsid);

			// get the dimension of the DEM and the resolution
			height = demWR.getHeight();
			width = demWR.getWidth();
			dx = CoverageUtilities.getRegionParamsFromGridCoverage(inDem).get(CoverageUtilities.XRES);

			// compute the vector normal to a grid cell surface.
			normalWR = normalVector(demWR, dx);
		}

		// create the set of the coordinate of the station, so we can
		// iterate over the set
		Iterator<Integer> idIterator = stationCoordinates.keySet().iterator();

		// trasform the list of idStation into an array
		idStations = stationCoordinates.keySet().toArray();

		// iterate over the list of the stations to detect their position in the
		// map and their latitude
		// iterate over the list of the stations
		for (int i = 0; i < idStations.length; i++) {

			// compute the coordinate of the station from the linked hashMap
			Coordinate coordinate = (Coordinate) stationCoordinates.get(idIterator.next());

			// define the position, according to the CRS, of the station in the map
			DirectPosition point = new DirectPosition2D(sourceCRS, coordinate.x, coordinate.y);

			// trasform the position in two the indices of row and column
			DirectPosition gridPoint = transf.transform(point, null);

			// add the indices to a list
			columnStation.add((int) gridPoint.getCoordinate()[0]);
			rowStation.add((int) gridPoint.getCoordinate()[1]);

			// reproject the map in WGS84 and compute the latitude
			Point[] idPoint = getPoint(coordinate, sourceCRS, targetCRS);
			latitudeStation.add(Math.toRadians(idPoint[0].getY()));
			longitudeStation.add(Math.toRadians(idPoint[0].getX()));

			// read the input data for the given station
			temperature = defaultTemp;
			if (inTemperatureValues != null)
				temperature = inTemperatureValues.get(idStations[i])[0];

			humidity = pRH;
			if (inHumidityValues != null)
				humidity = inHumidityValues.get(idStations[i])[0];

			// calculating the sun vector
			double sunVector[] = calcSunVector(latitudeStation.get(i),
					getHourAngle(date, latitudeStation.get(i), longitudeStation.get(i)));

			// calculate the inverse of the sun vector
			double[] inverseSunVector = calcInverseSunVector(sunVector);

			// calculate the normal of the sun vector according to Corripio (2003)
			double[] normalSunVector = calcNormalSunVector(sunVector);

			// E0 is the correction factor related to Earth’s orbit eccentricity computed
			// according to Spencer (1971):
			double E0 = computeE0(date);

			// evaluate the shadow map
			WritableRaster shadowWR = calculateFactor(height, width, sunVector, inverseSunVector, normalSunVector,
					demWR, dx);

			// calculate the direct radiation, during the daylight
			double direct = (hour > (sunrise) && hour < (sunset))
					? calcDirectRadiation(columnStation.get(i), rowStation.get(i), demWR, shadowWR, sunVector, normalWR,
							E0, temperature, humidity)
					: 0;

			// calculate the diffuse radiation, during the daylight
			double diffuse = (hour > (sunrise) && hour < (sunset))
					? calcDiffuseRadiation(sunVector, E0, columnStation.get(i), rowStation.get(i))
					: 0;

			// calculate the radiation at the top of the atmosphere, during the daylight
			double topATM = (hour > (sunrise) && hour < (sunset)) ? calcTopAtmosphere(E0, sunVector[2]) : 0;

			// calculate the radiation at the top of the atmosphere, during the daylight
			double total = (hour > (sunrise) && hour < (sunset)) ? direct + diffuse : 0;

			storeResult_series((Integer) idStations[i], direct, diffuse, topATM, total);

		}

		// upgrade the step for the date
		step++;
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
		return 1.00011 + 0.034221 * Math.cos(k) + 0.00128 * Math.sin(k) + 0.000719 * Math.cos(2 * k)
				+ 0.000077 * Math.sin(2 * k);
	}

	/**
	 * Gets the coordinate given the shp file and the field name in the shape with
	 * the coordinate of the station.
	 *
	 * @param collection is the shp file with the stations
	 * @param idField    is the name of the field with the id of the stations
	 * @return the coordinate of each station
	 * @throws Exception the exception in a linked hash map
	 */
	private LinkedHashMap<Integer, Coordinate> getCoordinate(SimpleFeatureCollection collection, String idField)
			throws Exception {
		LinkedHashMap<Integer, Coordinate> id2CoordinatesMap = new LinkedHashMap<Integer, Coordinate>();
		FeatureIterator<SimpleFeature> iterator = collection.features();
		Coordinate coordinate = null;
		try {
			while (iterator.hasNext()) {
				SimpleFeature feature = iterator.next();
				int stationNumber = ((Number) feature.getAttribute(idField)).intValue();
				coordinate = ((Geometry) feature.getDefaultGeometry()).getCentroid().getCoordinate();
				id2CoordinatesMap.put(stationNumber, coordinate);
			}
		} finally {
			iterator.close();
		}

		return id2CoordinatesMap;
	}

	/**
	 * Gets the point coordinates (row and column) of the station, after its
	 * reprojection in WGS84.
	 *
	 * @param coordinate is the coordinate of the point in the original reference
	 *                   system
	 * @param sourceCRS  the original reference system
	 * @param targetCRS  is the WGS84 system
	 * @return the point vector with the x and y values of its position
	 * @throws Exception
	 */
	private Point[] getPoint(Coordinate coordinate, CoordinateReferenceSystem sourceCRS,
			CoordinateReferenceSystem targetCRS) throws Exception {
		Point[] point = new Point[] { GeometryUtilities.gf().createPoint(coordinate) };
		CrsUtilities.reproject(sourceCRS, targetCRS, point);
		return point;
	}

	/**
	 * getHourAngle is the value of the hour angle at a given time and latitude
	 * (Corripio (2003))
	 * 
	 *
	 * @param date     is the current date
	 * @param latitude is the latitude of the station
	 * @return the double value of the hour angle
	 */

	public double radians(double degree) {
		double radian = degree * (Math.PI / 180.0);
		return (radian);
	}

	public double degrees(double radian) {
		double degree = radian * (180.0 / Math.PI);
		return (degree);
	}

	private double getHourAngle(DateTime date, double latitude, double longitude) {
		int day = date.getDayOfYear();

		double jd = ((dateToJulian(date.toCalendar(null))));
		hour = (doHourly == true) ? (double) date.getMillisOfDay() / (1000 * 60 * 60) : 12.5;

		hour = ((jd - Math.floor(jd)) * 24 + 12) % 24;

		double jdc = (jd - 2451545.0) / 36525.0;
		double sec = 21.448 - jdc * (46.8150 + jdc * (0.00059 - jdc * (0.001813)));
		double e0 = 23.0 + (26.0 + (sec / 60.0)) / 60.0;
		double ecc = 0.016708634 - jdc * (0.000042037 + 0.0000001267 * jdc);
		double oblcorr = e0 + 0.00256 * Math.cos(radians(125.04 - 1934.136 * jdc));
		double y = Math.pow((Math.tan(radians(oblcorr) / 2)), 2);
		double l0 = 280.46646 + jdc * (36000.76983 + jdc * (0.0003032));
		l0 = (l0 - 360 * ((int) l0 / (int) 360)) % 360;
		double rl0 = radians(l0);
		double gmas = 357.52911 + jdc * (35999.05029 - 0.0001537 * jdc);
		gmas = radians(gmas);
		double EqTime = y * Math.sin(2 * rl0) - 2.0 * ecc * Math.sin(gmas)
				+ 4.0 * ecc * y * Math.sin(gmas) * Math.cos(2 * rl0) - 0.5 * y * y * Math.sin(4 * rl0)
				- 1.25 * ecc * ecc * Math.sin(2 * gmas);

		double timezone = 0;
		double eqtime = degrees(EqTime) * 4;
		double stndmeridian = timezone * 15;
		double deltalontime = degrees(longitude) - stndmeridian;
		deltalontime = deltalontime * 24.0 / 360.0;
		double omegar = Math.PI * (((hour + deltalontime + eqtime / 60) / 12.0) - 1.0);

		// Evaluate the declination of the sun.
		delta = radians(declination(jd));
		// ss is the absolute value of the hour angle at sunrise or sunset

		double tanlatdel = -Math.tan((latitude)) * Math.tan((delta));
		if (tanlatdel > 1)
			tanlatdel = 1;
		double omega = Math.acos(tanlatdel);
		double daylen = (2 * omega) / (2 * Math.PI / 24);
		sunrise = 12 * (1 - omega / Math.PI) - deltalontime - EqTime / 60;
		sunset = 12 * (1 + omega / Math.PI) - deltalontime - EqTime / 60;
		if (omega == 0)
			sunrise = -9999;
		if (omega == 0)
			sunset = -9999;

		double hourangle = omegar;
		return hourangle;

	}

	public double dateToJulian(Calendar date) {
		int year = date.get(Calendar.YEAR);
		int month = date.get(Calendar.MONTH) + 1;
		int day = date.get(Calendar.DAY_OF_MONTH);
		int hour = date.get(Calendar.HOUR_OF_DAY);
		int minute = date.get(Calendar.MINUTE);
		int second = date.get(Calendar.SECOND);

		double extra = (100.0 * year) + month - 190002.5;
		return (367.0 * year) - (Math.floor(7.0 * (year + Math.floor((month + 9.0) / 12.0)) / 4.0))
				+ Math.floor((275.0 * month) / 9.0) + day + ((hour + ((minute + (second / 60.0)) / 60.0)) / 24.0)
				+ 1721013.5 - ((0.5 * extra) / Math.abs(extra)) + 0.5;
	}

	public double declination(double jd) {

		double T = (jd - 2451545) / 36525.0;
		double epsilon = (23 + 26 / 60.0 + 21.448 / 3600.0) - (46.8150 / 3600.0) * T - (0.00059 / 3600.0) * T * T
				+ (0.001813 / 3600.0) * T * T * T;
		double L0 = 280.46645 + 36000.76983 * T + 0.0003032 * T * T;
		L0 = (L0 - 360 * ((int) L0 / (int) 360)) % 360;
		double M = 357.52910 + 35999.05030 * T - 0.0001559 * T * T - 0.00000048 * T * T * T;
		double e = 0.016708617 - 0.000042037 * T - 0.0000001236 * T * T;
		double C = (1.914600 - 0.004817 * T - 0.000014 * T * T) * Math.sin(radians(M))
				+ (0.019993 - 0.000101 * T) * Math.sin(2 * radians(M)) + 0.000290 * Math.sin(3 * radians(M));
		double Theta = L0 + C;
		double v = M + C;
		double Omega = 125.04452 - 1934.136261 * T + 0.0020708 * T * T + (T * T * T) / 450000;
		double lambda = Theta - 0.00569 - 0.00478 * Math.sin(radians(Omega));
		double delta = Math.asin(Math.sin(radians(epsilon)) * Math.sin(radians(lambda)));
		return (degrees(delta));

	}

	/**
	 * calcSunVector compute the vector vector in the direction of the Sun (Corripio
	 * (2003))
	 *
	 * @param latitude is the latitude of the station
	 * @param the      hour angle
	 * @return the sun vector
	 */
	protected double[] calcSunVector(double latitude, double hourAngle) {
		double sunVector[] = new double[3];
		sunVector[0] = -Math.sin(hourAngle) * Math.cos(delta);
		sunVector[1] = Math.sin((latitude)) * Math.cos(hourAngle) * Math.cos(delta)
				- Math.cos((latitude)) * Math.sin(delta);
		sunVector[2] = Math.cos((latitude)) * Math.cos(hourAngle) * Math.cos(delta)
				+ Math.sin((latitude)) * Math.sin(delta);
		return sunVector;

	}

	/**
	 * Maps reader transform the GrifCoverage2D in to the writable raster, replace
	 * the -9999.0 value with no value.
	 *
	 * @param inValues: the input map values
	 * @return the writable raster of the given map
	 */
	private WritableRaster mapsTransform(GridCoverage2D inValues) {
		RenderedImage inValuesRenderedImage = inValues.getRenderedImage();
		WritableRaster inValuesWR = CoverageUtilities.replaceNovalue(inValuesRenderedImage, -9999.0);
		inValuesRenderedImage = null;
		return inValuesWR;
	}

	/**
	 * normalVector compute the vector normal to a grid cell surface, according to
	 * Corripio (2003)
	 *
	 * @param demWR is the Writable raster of DEM
	 * @param res   is the resolution of the DEM
	 * @return the normal vector for each cell
	 */
	protected WritableRaster normalVector(WritableRaster demWR, double res) {

		int minX = demWR.getMinX();
		int minY = demWR.getMinY();
		int rows = demWR.getHeight();
		int cols = demWR.getWidth();

		RandomIter pitIter = RandomIterFactory.create(demWR, null);
		/*
		 * Initialize the image of the normal vector in the central point of the cells,
		 * which have 3 components (X;Y;Z), so the Image have 3 bands..
		 */
		SampleModel sm = RasterFactory.createBandedSampleModel(5, cols, rows, 3);
		WritableRaster tmpNormalVectorWR = CoverageUtilities.createDoubleWritableRaster(cols, rows, null, sm, 0.0);
		WritableRandomIter tmpNormalIter = RandomIterFactory.createWritable(tmpNormalVectorWR, null);
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
				double den = Math.sqrt(firstComponent * firstComponent + secondComponent * secondComponent
						+ thirthComponent * thirthComponent);
				tmpNormalIter.setPixel(i, j,
						new double[] { firstComponent / den, secondComponent / den, thirthComponent / den });

			}
		}
		pitIter.done();

		return tmpNormalVectorWR;

	}

	/**
	 * calcDirectRadiation calculates the direct radiation according to Corripio
	 * (2002)
	 *
	 * @param i          is the the index of the column of the station in the DEM
	 * @param j          is the the index of the row of the station in the DEM
	 * @param demWR      is the DEM
	 * @param shadowWR   is the map with the shadow
	 * @param sunVector  is the sun vector
	 * @param normalWR   is the the raster with the vector normal to a pixel
	 * @param E0         is the correction of the eccentricity
	 * @param temperture is the air temperature
	 * @param himidity   is the relative humidity
	 * @return the double value of the direct radiation
	 */

	private double calcDirectRadiation(int i, int j, WritableRaster demWR, WritableRaster shadowWR, double[] sunVector,
			WritableRaster normalWR, double E0, double temperature, double humidity) throws IOException {

		// zenith angle
		double zenith = Math.acos(sunVector[2]);

		// mr [–] relative optical air mass:
		double mr = 1.0 / (sunVector[2] + 0.15 * Math.pow((93.885 - (zenith * (180 / (2 * Math.PI)))), (-1.253)));
		mr = 1.0 / (Math.cos(zenith) + 0.15 * ((Math.pow(93.885 - degrees(zenith), -1.253))));

		// altitude of the station
		double z = demWR.getSampleDouble(i, j, 0);

		// local atmospheric pressure
		double pressure = ATM * Math.exp(-0.0001184 * z);

		// relative air mass
		ma = mr * pressure / ATM;

		// The transmittance functions for Rayleigh scattering
		tau_r = Math.exp((-.09030 * Math.pow(ma, 0.84)) * (1.0 + ma - Math.pow(ma, 1.01)));

		// transform the temperature in Kelvin
		temperature = temperature + 273.0;

		// evaluate the saturated valor pressure
		double saturatedVaporPressure = Math.exp(26.23 - 5416.0 / temperature);

		// the precipitable water in cm calculated according to Prata (1996)
		double w = 0.465 * (humidity) * saturatedVaporPressure / temperature;

		// Transmittance by water vapour
		tau_w = 1.0 - 2.4959 * w * mr / (Math.pow(1.0 + 79.034 * w * mr, 0.6828) + 6.385 * w * mr);

		// The transmittance by ozone
		tau_o = 1.0 - ((0.1611 * pCmO3 * mr * Math.pow(1.0 + 139.48 * pCmO3 * mr, -0.3035))
				- (0.002715 * pCmO3 * mr / (1.0 + 0.044 * pCmO3 * mr + 0.0003 * Math.pow(pCmO3 * mr, 2))));

		// Transmittance by uniformly mixed gases
		tau_g = Math.exp(-0.0127 * Math.pow(ma, 0.26));

		// The transmittance by aerosols
		tau_a = Math.pow((0.97 - 1.265 * Math.pow(pVisibility, (-0.66))), Math.pow(ma, 0.9));

		// correction factor [m] for increased trasmittance with elevation z[m]
		// according to Corripio (2002)
		double beta_s = (z <= 3000) ? 2.2 * Math.pow(10, -5) * z : 2.2 * Math.pow(10, -5) * 3000;

		// cosin: incidence angle according to Corripio (2003)
		double cos_inc = scalarProduct(sunVector, normalWR.getPixel(i, j, new double[3]));
		if (cos_inc < 0)
			cos_inc = 0;

		// Direct radiation under cloudless sky incident on arbitrary tilted
		// surfaces (by inclusion of cosinc)
		In = 0.9571 * SOLARCTE * E0 * (tau_r * tau_o * tau_g * tau_w * tau_a + beta_s);

		double S_incident = In * cos_inc * shadowWR.getSampleDouble(i, j, 0);

		S_incident = (S_incident > 3000) ? 0 : S_incident;
		S_incident = (S_incident <= 0) ? 0 : S_incident;

		return S_incident;
	}

	/**
	 * calcDiffuseRadiation calculates the diffuse radiation according to Corripio
	 * (2002)
	 *
	 * @param i         is the the index of the column of the station in the DEM
	 * @param j         is the the index of the row of the station in the DEM
	 * @param sunVector is the sun vector
	 * @param E0        is the correction of the eccentricity
	 * @return the double value of the diffuse radiation
	 */
	private double calcDiffuseRadiation(double[] sunVector, double E0, int i, int j) {

		// single-scattering albedo fraction of incident energy scattered to total
		// attenuation by aerosol
		double ssctalb = 0.9;// single scattering albedo (aerosols)(Iqbal, 1983)
		double zenith = Math.acos(sunVector[2]);

		double tauaa = 1.0 - (1.0 - ssctalb) * (1.0 - ma + Math.pow(ma, 1.06)) * (1.0 - tau_a);
		double I_dr = 0.79 * E0 * SOLARCTE * Math.cos(zenith) * tau_o * tau_g * tau_w * tauaa * 0.5 * (1.0 - tau_r)
				/ (1.0 - ma + Math.pow(ma, 1.02));
		double FC = 0.74;

		double tauas = (tau_a) / tauaa;
		double I_da = 0.79 * E0 * SOLARCTE * Math.cos(zenith) * tau_o * tau_g * tau_w * tauaa * FC * (1.0 - tauas)
				/ (1.0 - ma + Math.pow(ma, 1.02));

		// The atmospheric albedo is computed as
		double alpha_a = 0.0685 + (1.0 - FC) * (1.0 - (tau_a / tauaa));

		// the diffuse irradiance from multiple reflection between the earth and the
		// atmosphere
		double I_dm = (In * sunVector[2] + I_dr + I_da) * alpha_a * pAlphag / (1.0 - pAlphag * alpha_a);

		double diffuse = (I_dr + I_da + I_dm) * skyviewfactorWR.getSampleDouble(i, j, 0);

		diffuse = (diffuse > 3000) ? 0 : diffuse;
		diffuse = (diffuse < 0) ? 0 : diffuse;

		return diffuse;

	}

	/**
	 * calcTopAtmosphere calculates the radiation at the top of the atmosphere
	 * according to Corripio (2002)
	 *
	 * @param sunVector is the sun vector
	 * @param E0        is the correction of the eccentricity
	 * @return the double value of the top atmosphere radiation
	 */
	private double calcTopAtmosphere(double E0, double sunVector) {
		double topATM = E0 * SOLARCTE * Math.cos(Math.acos(sunVector));
		return topATM = (topATM < 0) ? 0 : topATM;
	}

	/**
	 * Store result_series stores the results in the hashMaps .
	 *
	 * @param ID      is the id of the station
	 * @param direct  is the direct radiation
	 * @param diffuse is the diffuse radiation
	 * @param topATM  is the radiation at the top of the atmosphere
	 * @throws SchemaException
	 */
	private void storeResult_series(int ID, double direct, double diffuse, double topATM, double total)
			throws SchemaException {
		outHMdirect.put(ID, new double[] { direct });

		outHMdiffuse.put(ID, new double[] { diffuse });

		outHMtopatm.put(ID, new double[] { topATM });

		outHMtotal.put(ID, new double[] { total });
	}

}