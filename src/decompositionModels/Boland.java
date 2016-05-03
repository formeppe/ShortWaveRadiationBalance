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
package decompositionModels;


// TODO: Auto-generated Javadoc
/**
 * The Class Boland implements the Boland [2001] model .
 */
public class Boland implements Model{


	/** The cleraness index. */
	double cleranessIndex;
	
	/** The kd. */
	double kd;


	/**
	 * Instantiates a new boland model.
	 *
	 * @param cleranessIndex the cleraness index
	 */
	public Boland(double cleranessIndex){

		this.cleranessIndex=cleranessIndex;
	}




	/* (non-Javadoc)
	 * @see decompositionModels.Model#kdValues()
	 */
	@Override
	public double kdValues() {
		
		kd = 1 / (1 + Math.exp(7.997 * (cleranessIndex - 0.586)));
		return kd;
	}


}
