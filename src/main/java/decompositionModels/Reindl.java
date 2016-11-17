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
 * The Class Reindl implemtes the Reindl [1990] model.
 */
public class Reindl implements Model{


	/** The cleraness index. */
	double cleranessIndex;
	
	/** The kd. */
	double kd;


	/**
	 * Instantiates a new reindl model.
	 *
	 * @param cleranessIndex the cleraness index
	 */
	public Reindl(double cleranessIndex){

		this.cleranessIndex=cleranessIndex;
	}




	/* (non-Javadoc)
	 * @see decompositionModels.Model#kdValues()
	 */
	@Override
	public double kdValues() {
		
		if (cleranessIndex < 0.3) {
			kd = 1.02 - 0.248 * cleranessIndex;
		} else if (cleranessIndex < 0.78) {
			kd = 1.45 - 1.67 * cleranessIndex;
		} else {
			kd = 0.147;
		}
		return kd;
	}


}
