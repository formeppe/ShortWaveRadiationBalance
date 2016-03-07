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

import org.geotools.coverage.grid.GridCoverage2D;
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
@Documentation("Insolation.html")
@Author(name = "Giuseppe Formettam Daniele Andreis and Riccardo Rigon", contact = "http://www.ing.unitn.it/dica/hp/?user=rigon")
@Keywords("Hydrology, Radiation, SkyviewFactor, Hillshade")
@Bibliography("Corripio, J. G.: 2003,"
		+ " Vectorial algebra algorithms for calculating terrain parameters"
		+ "from DEMs and the position of the sun for solar radiation modelling in mountainous terrain"
		+ ", International Journal of Geographical Information Science 17(1), 1–23. and"
		+ "Iqbal, M., 1983. An Introduction to solar radiation. In: , Academic Press, New York")
@Label(JGTConstants.HYDROGEOMORPHOLOGY)
@Name("shortradbal")
@Status(Status.CERTIFIED)
@License("General Public License Version 3 (GPLv3)")
public class ShortwaveRadiationBalanceRasterCase extends JGTModel {
	@Description("The map of the elevation.")
	@In
	public GridCoverage2D inElev = null;

	@Description("The map of the skyview factor.")
	@In
	public GridCoverage2D inskyview = null;

	@Description("The vector of the measurement points, containing the position of the stations.")
	@In
	public SimpleFeatureCollection inStations = null;

	@Description("The field of the vector of stations, defining the id.")
	@In
	public String fStationsid = null;

	@Description("The field of the vector of stations, defining the albedo.")
	@In
	public String fStationsAlb = null;

	@Description("The first day of the simulation.")
	@In
	public String tStartDate = null;

	@Description("The last day of the simulation.")
	@In
	public String tEndDate = null;
	@Description("The time step in minutes of the measurement.")
	@In
	public int inTimestep;

	@Description("The path to diffuse output.")
	@In
	public String pOutPathdiffuse;

	@Description("The path to tau output.")
	@In
	public String pOutPathttau;

	@Description("The path to direct output.")
	@In
	public String pOutPathdiretta;

	@Description("The path to top atm output.")
	@In
	public String pOutPathtopatm;

	@Description("Raster Mode=true; Vector Mode=F.")
	@In
	public boolean doRaster = true;

	@Description("The calculation time.")
	@In
	public int cumulationTime = 0;

	@Description("Path to the temperature file in input.")
	@In
	public String inPathtemp = null;

	@Description("Ozone layer thickness in cm.")
	@In
	public double pCmO3 = 0.4;

	@Description("Default relative humidity value.")
	@In
	public double pRH = 0.7;

	@Description(" For aerosol attenuation (5 < vis < 180 Km) [km].")
	@In
	public double pVisibility = 80;

	@Description("Path to the humidity file in input.")
	@In
	public String inPathhumidity = null;

	@Description("The progress monitor.")
	@In
	public IJGTProgressMonitor pm = new LogProgressMonitor();

	@Description("The map of total insolation in case workWithRaster=T.")
	@Out
	public GridCoverage2D outIns;

	@Description("The soil albedo.")
	@Out
	public double pAlphag = 0.9;

	// private static final double pCmO3 = 0.6;

	private static final double pLapse = -.0065;

	// ok private static final double pVisibility = 60;

	/**
	 * The solar constant.
	 */
	// ok private static final double SOLARCTE = 1360.0;
	private static final double SOLARCTE = 1370.0;
	/**
	 * The atmosphere pressure.
	 */
	private static final double ATM = 1013.25;

	private double lambda;
	private HashMap<Integer, double[]> outHMdirect;
	private HashMap<Integer, double[]> outHMdiffuse;
	private HashMap<Integer, double[]> outHMtopatm;

	private HashMap<String, Double> attribute;
	private double delta;
	private int height = 0;
	private int width = 0;
	private double omega;
	private double vetgio[];
	private int contaore = 1;
	private WritableRaster resultstaoWR = null;
	private WritableRaster resultinsWR = null;
	private WritableRaster resultdiffWR = null;
	private WritableRaster pitWR = null;
	private int contastampe = 0;
	private double[] xStation;
	private double[] yStation;
	private double[] zStation;
	private double[] albedoStation;

	private int[] idStation;
	private int[] colnumvetVect;
	private int[] rownumvetVect;
	private HashMap<Integer, double[]> tempvalues;
	private HashMap<Integer, double[]> umivalues;
	private WritableRaster skyviewfactorWR = null;

	@Execute
	public void process() throws Exception { // transform the

		// extract some attributes of the map
		attribute = CoverageUtilities.getRegionParamsFromGridCoverage(inElev);
		double dx = attribute.get(CoverageUtilities.XRES);

		/*
		 * The models use only one value of the latitude. So I have decided to
		 * set it to the center of the raster. Extract the CRS of the
		 * GridCoverage and transform the value of a WGS84 latitude.
		 */
		CoordinateReferenceSystem sourceCRS = inElev
				.getCoordinateReferenceSystem2D();
		CoordinateReferenceSystem targetCRS = DefaultGeographicCRS.WGS84;

		DateTimeFormatter formatter = DateTimeFormat.forPattern(
				"yyyy-MM-dd HH:mm").withZone(DateTimeZone.UTC);
		DateTime startcurrentDatetime = formatter.parseDateTime(tStartDate);
		DateTime startcurrentDatetimefittizio = startcurrentDatetime
				.plus(30 * 60 * 1000);

		DateTime endcurrentDatetime = formatter.parseDateTime(tEndDate);

		long diff = (endcurrentDatetime.getMillis() - startcurrentDatetime
				.getMillis())
				/ (1000 * 60 * 60);
		DateTime array[] = new DateTime[(int) diff];
		for (int i = 0; i < array.length; i++) {
			array[i] = startcurrentDatetimefittizio.plusHours(i);
		}
		if (doRaster == false) {
			List<Double> xStationList = new ArrayList<Double>();
			List<Double> yStationList = new ArrayList<Double>();
			List<Double> zStationList = new ArrayList<Double>();
			List<Double> iddList = new ArrayList<Double>();
			List<Double> albedoList = new ArrayList<Double>();

			/*
			 * counter for the number of station with measured value !=0.
			 */
			/*
			 * Store the station coordinates and measured data in the array.
			 */
			FeatureIterator<SimpleFeature> stationsIter = inStations.features();
			try {
				while (stationsIter.hasNext()) {
					SimpleFeature feature = stationsIter.next();
					double zzz = 0;
					int id = ((Number) feature.getAttribute(fStationsid))
							.intValue();
					if (fStationsAlb != null) {
						double albedo = ((Number) feature
								.getAttribute(fStationsAlb)).doubleValue();
						albedoList.add(albedo);
					}
					Coordinate coordinate = ((Geometry) feature
							.getDefaultGeometry()).getCentroid()
							.getCoordinate();
					xStationList.add(coordinate.x);
					yStationList.add(coordinate.y);
					zStationList.add(zzz);

					iddList.add((double) id);

				}
			} finally {
				stationsIter.close();
			}

			int nStaz = xStationList.size();
			/*
			 * The coordinates of the station points plus in last position a
			 * place for the coordinate of the point to interpolate.
			 */
			xStation = new double[nStaz];
			yStation = new double[nStaz];
			zStation = new double[nStaz];
			idStation = new int[nStaz];
			colnumvetVect = new int[nStaz];
			rownumvetVect = new int[nStaz];
			albedoStation = new double[nStaz];
			if (nStaz != 0) {
				for (int i = 0; i < nStaz; i++) {
					double xTmp = xStationList.get(i);
					double yTmp = yStationList.get(i);
					double zTmp = zStationList.get(i);
					int idTmp = iddList.get(i).intValue();
					xStation[i] = xTmp;
					yStation[i] = yTmp;
					zStation[i] = zTmp;
					idStation[i] = idTmp;
					if (albedoList.size() != 0l) {
						albedoStation[i] = albedoList.get(i);
					} else {
						if (!isNovalue(pAlphag)) {
							albedoStation[i] = pAlphag;
						} else {
							albedoStation[i] = 0.5;
						}
					}
				}
			}

			MathTransform transf = inElev.getGridGeometry().getCRSToGrid2D();
			for (int i = 0; i < xStation.length; i++) {

				DirectPosition point = new DirectPosition2D(sourceCRS,
						xStation[i], yStation[i]);
				DirectPosition gridPoint = transf.transform(point, null);

				colnumvetVect[i] = (int) gridPoint.getCoordinate()[0];
				rownumvetVect[i] = (int) gridPoint.getCoordinate()[1];
				System.out.println(idStation[i] + "  "
						+ gridPoint.getCoordinate()[0] + " "
						+ gridPoint.getCoordinate()[1]);
			}
		}

		double srcPts[] = new double[] { attribute.get(CoverageUtilities.EAST),
				attribute.get(CoverageUtilities.SOUTH) };

		Coordinate source = new Coordinate(srcPts[0], srcPts[1]);
		Point[] so = new Point[] { GeometryUtilities.gf().createPoint(source) };
		CrsUtilities.reproject(sourceCRS, targetCRS, so);
		// the latitude value
		lambda = Math.toRadians(so[0].getY());
		CoverageUtilities.getRegionParamsFromGridCoverage(inElev);
		
		
		RenderedImage pitTmpRI = inElev.getRenderedImage();
		width = pitTmpRI.getWidth();
		height = pitTmpRI.getHeight();
		pitWR = CoverageUtilities.replaceNovalue(pitTmpRI, -9999.0);
		pitTmpRI = null;

		WritableRaster insolationWR = CoverageUtilities
				.createDoubleWritableRaster(width, height, null, pitWR
						.getSampleModel(), 0.0);

		WritableRaster staoWR = CoverageUtilities.createDoubleWritableRaster(
				width, height, null, pitWR.getSampleModel(), 0.0);

		WritableRaster diffuseWR = CoverageUtilities
				.createDoubleWritableRaster(width, height, null, pitWR
						.getSampleModel(), 0.0);
		WritableRaster gradientWR = normalVector(pitWR, dx);

		//
		vetgio = new double[array.length];

		if (doRaster) {
			resultstaoWR = CoverageUtilities.createDoubleWritableRaster(width,
					height, null, pitWR.getSampleModel(), 0.0);
			resultinsWR = CoverageUtilities.createDoubleWritableRaster(width,
					height, null, pitWR.getSampleModel(), 0.0);
			resultdiffWR = CoverageUtilities.createDoubleWritableRaster(width,
					height, null, pitWR.getSampleModel(), 0.0);

			for (int i = 0; i < array.length; i++) {
				DateTime currentime = array[i];
				calcInsolation(lambda, pitWR, gradientWR, insolationWR, staoWR,
						diffuseWR, dx, currentime);

				// pm.worked(i - startDay);
			}
		}
		if (doRaster == false) {

			OmsTimeSeriesIteratorReader reader_rain = new OmsTimeSeriesIteratorReader();
			if (!((inPathtemp == null))) {
				reader_rain.file = inPathtemp;
				reader_rain.idfield = "ID";
				reader_rain.tStart = tStartDate;
				reader_rain.tEnd = tEndDate;
				reader_rain.fileNovalue = "-9999";
				reader_rain.tTimestep = inTimestep;
			}
			OmsTimeSeriesIteratorReader reader_umi = new OmsTimeSeriesIteratorReader();

			if (!(inPathhumidity == null)) {

				reader_umi.file = inPathhumidity;
				reader_umi.idfield = "ID";
				reader_umi.tStart = tStartDate;
				reader_umi.tEnd = tEndDate;
				reader_umi.fileNovalue = "-9999";
				reader_umi.tTimestep = inTimestep;

			}

			OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
			OmsTimeSeriesIteratorWriter writer2 = new OmsTimeSeriesIteratorWriter();
			OmsTimeSeriesIteratorWriter writer3 = new OmsTimeSeriesIteratorWriter();

			writer.file = pOutPathdiretta;
			writer2.file = pOutPathdiffuse;
			writer3.file = pOutPathtopatm;

			writer.tStart = tStartDate;
			writer2.tStart = tStartDate;
			writer3.tStart = tStartDate;
			writer2.tTimestep = inTimestep;
			writer.tTimestep = inTimestep;
			writer3.tTimestep = inTimestep;

			for (int i = 0; i < array.length; i++) {
				outHMdirect = new HashMap<Integer, double[]>();
				outHMdiffuse = new HashMap<Integer, double[]>();
				outHMtopatm = new HashMap<Integer, double[]>();
				System.out.println(" data=" + array[i]);
				DateTime currentime = array[i];
				if (!(inPathtemp == null)) {
					reader_rain.nextRecord();
					tempvalues = reader_rain.outData;
				}

				if (!(inPathhumidity == null)) {
					reader_umi.nextRecord();
					umivalues = reader_umi.outData;
				}

				calcInsolation(lambda, pitWR, gradientWR, insolationWR, staoWR,
						diffuseWR, dx, currentime);
				writer.inData = outHMdirect;
				writer.writeNextLine();
				writer2.inData = outHMdiffuse;
				writer2.writeNextLine();
				writer3.inData = outHMtopatm;
				writer3.writeNextLine();
				// pm.worked(i - startDay);
			}
			writer.close();
			writer2.close();
			writer3.close();
		}

	}

	/**
	 * Evaluate the radiation.
	 * 
	 * @param lambda
	 *            the latitude.
	 * @param demWR
	 *            the raster of elevation
	 * @param gradientWR
	 *            the raster of the gradient value of the dem.
	 * @param insolationWR
	 *            the wr where to store the result.
	 * @param the
	 *            day in the year.
	 * @throws Exception
	 * @paradx the resolutiono of the dem.
	 */
	private void calcInsolation(double lambda, WritableRaster demWR,
			WritableRaster gradientWR, WritableRaster insolationWR,
			WritableRaster staoWR, WritableRaster diffuseWR, double dx,
			DateTime time) throws Exception {

		int day = time.getDayOfYear();

		double dayangb = (360 / 365.25) * (day - 79.436);

		dayangb = Math.toRadians(dayangb);

		// Evaluate the declination of the sun.
		delta = getDeclination(dayangb);
		// Evaluate the radiation in this day.
		double ss = Math.acos(-Math.tan(delta) * Math.tan(lambda));
		double sunrise = 12 * (1.0 - ss / Math.PI);
		double sunset = 12 * (1.0 + ss / Math.PI);
		double hhh = (double) time.getMillisOfDay() / (1000 * 60 * 60);

		double tetaj = 2 * Math.PI * (day - 1.0) / 365.0;
		double ecccorr = 1.00011 + 0.034221 * Math.cos(tetaj) + 0.00128
				* Math.sin(tetaj) + 0.000719 * Math.cos(2 * tetaj) + 0.000077
				* Math.sin(2 * tetaj);

		double hourangle = (hhh / 12.0 - 1.0) * Math.PI;

		CoverageUtilities.getRegionParamsFromGridCoverage(inskyview);
		RenderedImage SkyTmpRI = inskyview.getRenderedImage();
		int width = SkyTmpRI.getWidth();
		int height = SkyTmpRI.getHeight();
		skyviewfactorWR = CoverageUtilities.renderedImage2WritableRaster(
				SkyTmpRI, true);

		if (doRaster) {

			for (int j = 0; j < height; j++) {
				for (int i = 0; i < width; i++) {

					if (hhh <= (sunrise)) {

						insolationWR.setSample(i, j, 0, 0);
						staoWR.setSample(i, j, 0, 0);
						diffuseWR.setSample(i, j, 0, 0);

					}
					if (hhh >= (sunset)) {
						// System.out.println(0);
						vetgio[contaore] = 0;
						insolationWR.setSample(i, j, 0, 0);
						staoWR.setSample(i, j, 0, 0);
						diffuseWR.setSample(i, j, 0, 0);

					}
					if (hhh > (sunrise) && hhh < (sunset)) {
						if (demWR.getSampleDouble(i, j, 0) != -9999.0) {
							if (hhh - sunrise < 0.01) {
								hhh = hhh + 0.1;
								hourangle = (hhh / 12.0 - 1.0) * Math.PI;
							}

							if (sunset - hhh < 0.01) {
								hhh = hhh - 0.1;
								hourangle = (hhh / 12.0 - 1.0) * Math.PI;
							}
							omega = hourangle;
							// calculating the vector related to the sun
							double sunVector[] = calcSunVector();
							double zenith = calcZenith(sunVector[2]);

							double[] inverseSunVector = calcInverseSunVector(sunVector);
							double[] normalSunVector = calcNormalSunVector(sunVector);

							height = demWR.getHeight();
							width = demWR.getWidth();
							WritableRaster sOmbraWR = calculateFactor(height,
									width, sunVector, inverseSunVector,
									normalSunVector, demWR, dx);

							double mr = 1 / (sunVector[2] + 0.15 * Math.pow(
									(93.885 - (zenith * 180 / 3.14)), (-1.253)));

							double temp = doubleNovalue;
							double umi = doubleNovalue;
							// evaluate the radiation.
							double[] aaa = calcRadiation(i, j, demWR, sOmbraWR,
									insolationWR, sunVector, gradientWR, mr,
									ecccorr, temp, umi);
							// vetgio[contaore] = aaa * 0.0864 / 24;
							// System.out.println(aaa);
							staoWR.setSample(i, j, 0, aaa[0] * 0.0864 / 24);
							insolationWR.setSample(i, j, 0,
									aaa[1] * 0.0864 / 24);
							diffuseWR.setSample(i, j, 0, aaa[2] * 0.0864 / 24);

						}

						else {
							staoWR.setSample(i, j, 0, -9999);
							insolationWR.setSample(i, j, 0, -9999);
							diffuseWR.setSample(i, j, 0, -9999);

						}
					}
				}
			}
			//
		} else {

			for (int contastaz = 0; contastaz < xStation.length; contastaz++) {
				System.out.println("STAZIONE N. ========> " + contastaz);
				int colnuber = colnumvetVect[contastaz];
				int rownumber = rownumvetVect[contastaz];
				int i = colnuber;
				int j = rownumber;
				int id = idStation[contastaz];
				if (hhh <= (sunrise)) {
					outHMdirect.put(id, new double[] { 0.0 });
					outHMdiffuse.put(id, new double[] { 0.0 });
					outHMtopatm.put(id, new double[] { 0.0 });
				}
				if (hhh >= (sunset)) {
					outHMdirect.put(id, new double[] { 0.0 });
					outHMdiffuse.put(id, new double[] { 0.0 });
					outHMtopatm.put(id, new double[] { 0.0 });
				}
				if (hhh > (sunrise) && hhh < (sunset)) {

					if (hhh - sunrise < 0.01) {
						hhh = hhh + 0.1;
						hourangle = (hhh / 12.0 - 1.0) * Math.PI;
					}

					if (sunset - hhh < 0.01) {
						hhh = hhh - 0.1;
						hourangle = (hhh / 12.0 - 1.0) * Math.PI;
					}
					omega = hourangle;

					// calculating the vector related to the sun
					double sunVector[] = calcSunVector();
					double zenith = calcZenith(sunVector[2]);

					double[] inverseSunVector = calcInverseSunVector(sunVector);
					double[] normalSunVector = calcNormalSunVector(sunVector);

					height = demWR.getHeight();
					width = demWR.getWidth();

					WritableRaster sOmbraWR = calculateFactor(height, width,
							sunVector, inverseSunVector, normalSunVector,
							demWR, dx);
					// //c'era e funzionava...lo tolgo perche il cos puo essere
					// negativo
					// if (sunVector[2] < 0) {
					// sunVector[2] = 0;
					// }
					// //c'era e funzionava...lo tolgo perche il cos puo essere
					// negativo
					
					// attenta qui perchè dovrebbe essere due pigreco 

					double mr = 1.0 / (sunVector[2] + 0.15 * Math.pow(
							(93.885 - (zenith * (180 / (2*3.14)))), (-1.253)));

					// evaluate the radiation.
					double temperatura = doubleNovalue;
					if (tempvalues != null) {
						temperatura = tempvalues.get(id)[0];
					}

					double umidita = doubleNovalue;
					if (umivalues != null) {
						umidita = umivalues.get(id)[0];
					}
					System.out.println("ID===" + id);
					double[] aaa = calcRadiation(i, j, demWR, sOmbraWR,
							insolationWR, sunVector, gradientWR, mr, ecccorr,
							temperatura, umidita);
					outHMdirect.put(id, new double[] { aaa[1] });
					outHMdiffuse.put(id, new double[] { aaa[2] });
					outHMtopatm.put(id, new double[] { aaa[0] });

				}
			}
		}

		//
		System.out.println("========contaora========= " + contaore);
		contaore += 1;
		if (doRaster) {
			printmap(insolationWR, staoWR, diffuseWR, cumulationTime, contaore);
		}
		// System.out.println("___________________");
	}

	public void printmap(WritableRaster rasterIsolation,
			WritableRaster rasterstao, WritableRaster rasterdiffuse,
			int cumulata, int ore) throws Exception {

		if (ore <= cumulata) {
			// System.out.println(ore+"  "+contaore);
			for (int j = 0; j < height; j++) {
				for (int i = 0; i < width; i++) {
					// evaluate the radiation.
					double valueins = rasterIsolation.getSampleDouble(i, j, 0)
							+ resultinsWR.getSampleDouble(i, j, 0);
					double valuestau = rasterstao.getSampleDouble(i, j, 0)
							+ resultstaoWR.getSampleDouble(i, j, 0);
					double valuediff = rasterdiffuse.getSampleDouble(i, j, 0)
							+ resultdiffWR.getSampleDouble(i, j, 0);
					// System.out.println(i+"  "+j+"  "+value);
					// System.out.println("somma=" + value);
					resultstaoWR.setSample(i, j, 0, valuestau);
					resultdiffWR.setSample(i, j, 0, valuediff);
					resultinsWR.setSample(i, j, 0, valueins);

				}
			}
		} else {

			for (int y = 2; y < height - 2; y++) {
				for (int x = 2; x < width - 2; x++) {
					if (pitWR.getSampleDouble(x, y, 0) == -9999.0) {
						resultdiffWR.setSample(x, y, 0, Double.NaN);
						resultinsWR.setSample(x, y, 0, Double.NaN);
						resultstaoWR.setSample(x, y, 0, Double.NaN);
					}

				}
			}
			contaore = 0;

			contastampe += 1;
			String resstau = "map_STAU" + contastampe;
			String percorsostau = "/Users/giuseppeformetta/Desktop/ValidationJAMILittleWashita/"
					+ resstau + ".asc";
			String resdiff = "map_DIFF" + contastampe;
			String percorsodiff = "/Users/giuseppeformetta/Desktop/ValidationJAMILittleWashita/"
					+ resdiff + ".asc";
			String resins = "map_INS" + contastampe;
			String percorsoins = "/Users/giuseppeformetta/Desktop/ValidationJAMILittleWashita/"
					+ resins + ".asc";

			GridCoverage2D coveragestau = CoverageUtilities.buildCoverage(
					"resultstao", resultstaoWR, attribute, inElev
							.getCoordinateReferenceSystem());
			GridCoverage2D coverageins = CoverageUtilities.buildCoverage(
					"resultins", resultinsWR, attribute, inElev
							.getCoordinateReferenceSystem());
			GridCoverage2D coveragediff = CoverageUtilities.buildCoverage(
					"resultdiff", resultdiffWR, attribute, inElev
							.getCoordinateReferenceSystem());

			OmsRasterWriter writer1 = new OmsRasterWriter();
			writer1.file = percorsostau;
			writer1.inRaster = coveragestau;
			writer1.process();

			OmsRasterWriter writer2 = new OmsRasterWriter();
			writer2.file = percorsoins;
			writer2.inRaster = coverageins;
			writer2.process();

			OmsRasterWriter writer3 = new OmsRasterWriter();
			writer3.file = percorsodiff;
			writer3.inRaster = coveragediff;
			writer3.process();

			resultdiffWR = null;
			resultdiffWR = CoverageUtilities.createDoubleWritableRaster(width,
					height, null, pitWR.getSampleModel(), 0.0);

			resultinsWR = null;
			resultinsWR = CoverageUtilities.createDoubleWritableRaster(width,
					height, null, pitWR.getSampleModel(), 0.0);

			resultstaoWR = null;
			resultstaoWR = CoverageUtilities.createDoubleWritableRaster(width,
					height, null, pitWR.getSampleModel(), 0.0);

		}

	}

	/*
	 * Evaluate the declination.
	 */
	private double getDeclination(double dayangb) {
		double delta = .3723 + 23.2567 * Math.sin(dayangb) - .758
				* Math.cos(dayangb) + .1149 * Math.sin(2 * dayangb) + .3656
				* Math.cos(2 * dayangb) - .1712 * Math.sin(3 * dayangb) + .0201
				* Math.cos(3 * dayangb);
		return Math.toRadians(delta);
	}

	/*
	 * evaluate several component of the radiation and then multiply by the
	 * sOmbra factor.
	 */
	private double[] calcRadiation(int i, int j, WritableRaster demWR,
			WritableRaster sOmbraWR, WritableRaster insolationWR,
			double[] sunVector, WritableRaster gradientWR, double mr,
			double eccentricita, double tem, double umi) throws IOException {

		double risultato[] = new double[3];
		if (demWR.getSampleDouble(i, j, 0) != -9999.0) {
			double diffuseshort = 0;
			double z = demWR.getSampleDouble(i, j, 0);
			double pressure = ATM * Math.exp(-0.0001184 * z);
			double ma = mr * pressure / ATM;
			//double psup0 = (288 - 0.0065 * z) / 288;
			//psup0 = Math.pow(psup0, 5.256);
			//double newMA = mr * psup0;
			//ma = newMA;
			double temp = tem + 273.0;
			double pRHH = umi;

			if (isNovalue(umi)) {
				pRHH = pRH;
			}

			if (isNovalue(tem)) {
				temp = 273.0 + pLapse * (z - 4000);
			}

			double vap_psat = Math.exp(26.23 - 5416.0 / temp);
			double wPrec = 0.493 * (pRHH / 100) * vap_psat / temp;

			// The transmittance functions for Rayleigh scattering
			double taur = Math.exp((-.09030 * Math.pow(ma, 0.84))
					* (1.0 + ma - Math.pow(ma, 1.01)));
			// ////////////////////////////////////////////////////

			// The transmittance by ozone // controllare bene questa formula lo
			// 0.0044
			double d = pCmO3 * mr;
			double tauo = 1.0 - ((0.1611 * d * Math.pow(1.0 + 139.48 * d,
					-0.3035)) - (0.002715 * d / (1.0 + 0.044 * d + 0.0003 * Math
					.pow(d, 2))));

			// ////////////////////////////////////////////////////

			// Transmittance by uniformly mixed gases
			double taug = Math.exp(-0.0127 * Math.pow(ma, 0.26));
			// ////////////////////////////////////////////////////

			// Transmittance -y water vapour
			double dd = wPrec * mr;
			double tauw = 1.0 - 2.4959 * (dd)
					/ (Math.pow(1.0 + 79.034 * dd, 0.6828) + 6.385 * (dd));

			// ////////////////////////////////////////////////////

			// The transmittance by aerosols
			double taua = Math.pow((0.97 - 1.265 * Math.pow(pVisibility,
					(-0.66))), Math.pow(ma, 0.9));
			// ////////////////////////////////////////////////////

			// daniele double In = 0.9751 * SOLARCTE * taur * tauo * taug * tauw
			// * taua;
			// mia vedi tesi corripio pag 22 formula 3.7 aggiungo beta(z)

			double betaz = 0;
			if (z <= 3000) {
				betaz = 2.2 * Math.pow(10, -5) * z;
			} else {
				betaz = 2.2 * Math.pow(10, -5) * 3000;
			}

			// ////////////////////////////////////////////////////

			// //calcolo secondo corripio dell cosinc
			double cosinc = scalarProduct(sunVector, gradientWR.getPixel(i, j,
					new double[3]));

			if (cosinc < 0) {
				cosinc = 0;
			}
			// ////////////////////////////////////////////////////

			// Direct radiation under cloudless sky incident on arbitrary tilted
			// surfaces (by inclusion of cosinc)

			double exoatmIrrad = eccentricita * SOLARCTE;
			double In = 0.9751 * exoatmIrrad
					* (taur * tauo * taug * tauw * taua + betaz);
			double rrr = In * cosinc * sOmbraWR.getSampleDouble(i, j, 0);

			// ////////////////////////////////////////////////////

			// ////////////////////////////////////////////////////
			// DIFFUSE COMPONENT
			// ////////////////////////////////////////////////////

			// Transmittance of direct radiation due to aerosol absorptance and
			// is given by

			double omegazero = 0.9;
			double tauaa = 1.0 - (1.0 - omegazero)
					* (1.0 - ma + Math.pow(ma, 1.06)) * (1 - taua);
			// ////////////////////////////////////////////////////

			// Components corresponding to the Rayleigh scattering after the
			// �rst
			double costeta = Math.min(sunVector[2], sunVector[2]);
			double Idr = 0.79 * exoatmIrrad * costeta * (1.0 - taur)
					* (tauo * taug * tauw * tauaa) * 0.5
					/ (1.0 - ma + Math.pow(ma, 1.02));

			// ////////////////////////////////////////////////////

			// The aerosol-scattered diffuse irradiance after the �rst pass
			// through the atmosphere is:
			// //per LW
			// double FC = 0.74;
			// //per LW
			double FC = 0.74;
			double Ida = 0.79 * exoatmIrrad * costeta
					* (tauo * taug * tauw * tauaa) * FC
					* (1.0 - (taua / tauaa)) / ((1 - ma + Math.pow(ma, 1.02)));

			// The atmospheric albedo is computed as

			double alphaa = 0.0685 + (1.0 - FC) * (1.0 - (taua / tauaa));

			// alphag e' l'albedo del suolo The value of is estimated from local
			// knowledge of the average
			// snow cover in a few km around the site of interest, as explained
			// by Greuell et al. (1997) and shown in

			double Idm = (In * sunVector[2] + Idr + Ida) * alphaa
					* pAlphag / (1.0 - pAlphag * alphaa);

			diffuseshort = (Idr + Ida + Idm)
					* skyviewfactorWR.getSampleDouble(i, j, 0);

			// System.out.println("albedo*= " + alphag * alphaa
			// + "  diffuseshort= " + diffuseshort);

			if (rrr > 3000) {
				rrr = doubleNovalue;
			}
			if (diffuseshort > 3000) {
				diffuseshort = doubleNovalue;
			}
			double zenithh = calcZenith(sunVector[2]);
			risultato[0] = eccentricita * SOLARCTE * Math.cos(zenithh);
			risultato[1] = rrr;
			risultato[2] = diffuseshort;

		}

		else {

			risultato[0] = -9999.0;
			risultato[1] = -9999.0;
			risultato[2] = -9999.0;

		}
		return risultato;

	}

	protected double[] calcSunVector() {
		double sunVector[] = new double[3];
		sunVector[0] = -Math.sin(omega) * Math.cos(delta);
		sunVector[1] = Math.sin(lambda) * Math.cos(omega) * Math.cos(delta)
				- Math.cos(lambda) * Math.sin(delta);
		sunVector[2] = Math.cos(lambda) * Math.cos(omega) * Math.cos(delta)
				+ Math.sin(lambda) * Math.sin(delta);

		return sunVector;

	}

	protected WritableRaster normalVector(WritableRaster pitWR, double res) {

		int minX = pitWR.getMinX();
		int minY = pitWR.getMinY();
		int rows = pitWR.getHeight();
		int cols = pitWR.getWidth();

		RandomIter pitIter = RandomIterFactory.create(pitWR, null);
		/*
		 * Initializa the Image of the normal vector in the central point of the
		 * cells, which have 3 components so the Image have 3 bands..
		 */
		SampleModel sm = RasterFactory
				.createBandedSampleModel(5, cols, rows, 3);
		WritableRaster tmpNormalVectorWR = CoverageUtilities
				.createDoubleWritableRaster(cols, rows, null, sm, 0.0);
		WritableRandomIter tmpNormalIter = RandomIterFactory.createWritable(
				tmpNormalVectorWR, null);
		/*
		 * apply the corripio's formula (is the formula (3) in the article)
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

	private double calcZenith(double sunVector2) {
		return Math.acos(sunVector2);
	}

}