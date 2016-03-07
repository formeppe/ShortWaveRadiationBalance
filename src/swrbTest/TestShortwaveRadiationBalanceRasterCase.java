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
package swrbTest;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureReader;

import swrbPointCase.ShortwaveRadiationBalance;
import swrbRasterCase.ShortwaveRadiationBalanceRasterCase;

import org.jgrasstools.hortonmachine.utils.HMTestCase;

/**
 * Test the {@link Insolation} module.
 * 
 * @author Giuseppe Formetta
 */
public class TestShortwaveRadiationBalanceRasterCase extends HMTestCase {

	private final static String START_DATE = "2002-01-01 07:00";
	private final static String END_DATE = "2002-01-07 00:00";

	public void testInsolation() throws Exception {

		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "/Users/marialaura/Desktop/ValidationJAMILittleWashita/stazioniGIUSTE.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		OmsRasterReader reader = new OmsRasterReader();
		reader.file = "/Users/marialaura/Desktop/LW/pit_LW.asc";
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D pit = reader.outRaster;

		OmsRasterReader readers = new OmsRasterReader();
		readers.file = "/Users/marialaura/Desktop/LW/sky.asc";
		readers.fileNovalue = -9999.0;
		readers.geodataNovalue = Double.NaN;
		readers.process();
		GridCoverage2D skyviewfactor = readers.outRaster;

		ShortwaveRadiationBalanceRasterCase insolation = new ShortwaveRadiationBalanceRasterCase();
		insolation.inStations = stationsFC;
		insolation.inElev = pit;
		insolation.inskyview = skyviewfactor;
		insolation.tStartDate = START_DATE;
		insolation.tEndDate = END_DATE;
		insolation.doRaster = false;
		insolation.fStationsid = "int_1";
		insolation.pOutPathdiffuse = "/Users/marialaura/Desktop/DIFFUSA.csv";
		insolation.pOutPathdiretta = "/Users/marialaura/Desktop/DIRETTA.csv";
		insolation.pOutPathtopatm = "/Users/marialaura/Desktop/TOPATMN.csv";
		insolation.inPathtemp = "/Users/marialaura/Desktop/ValidationJAMILittleWashita/temperature_orarie_2002_2008_NEWNEW.csv";
		insolation.inPathhumidity = "/Users/marialaura/Desktop/ValidationJAMILittleWashita/humidity_orarie_2002_2008_NEW.csv";
		insolation.inTimestep = 60;

		insolation.pm = pm;

		insolation.process();

	}

}
