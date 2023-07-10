import { Injectable } from '@angular/core';
import { Pool, RegionInfo, Region, PrimerPair, Primer } from '../types';
import resultData from '../../assets/results.json';
import regionData from '../../assets/regions.json';
import { groupBy } from '../utils';
import { CONFIG } from '../config';

@Injectable({
  providedIn: 'root'
})
export class ResultsService {
  private _error: boolean = false;
  private _results: Map<string, Region> = new Map<string, Region>();
  private _pools: Pool[] = (resultData as Pool[]);

  private _initResults(): void {
    // Load the results and regions data
    const regions: RegionInfo[] = regionData as RegionInfo[];

    // Concat each primer pair into a single list
    const listOfAllPrimers = this._pools.reduce((acc: PrimerPair[], pool: Pool) => {
      // Add pool_id to PrimerPair
      pool.primer_pairs.forEach((primerPair: PrimerPair) => primerPair.pool_id = pool.pool_id);

      // Add the primer pairs to the list
      return acc.concat(pool.primer_pairs);
    }, []);

    // Group the primer pairs by region name
    const results = new Map<string, Region>();
    groupBy<PrimerPair>(listOfAllPrimers, 'region_name')
      .forEach((primerPairs: PrimerPair[], key: string) => {
        // Find the region info for this region
        const regionInfo = regions.find((element: RegionInfo) => element.name === key);
        if (regionInfo === undefined) {
          throw new Error(CONFIG.MESSAGES.NO_REGION_INFO);
        }

        /*
        Group the primer pairs by pool_id, so that we can iterate over each region and pool_id
        to simply retrieve the primer pairs for that given region and pool_id. This way we can
        simply retrieve an amplicon by observing the forward_primer start and reverse_primer end
        */
        const primerPairsGroupedByPoolId =  groupBy<PrimerPair>(primerPairs, 'pool_id');
        primerPairsGroupedByPoolId.forEach((primerPairs: PrimerPair[]) => {
          // Sort the primer pairs by position using the forward primer
          primerPairs.sort((a: PrimerPair, b: PrimerPair) => a.forward_primer.position - b.forward_primer.position);

          // Asign random hex code color to the primer pair
          const validColors = CONFIG.COLORS;
          primerPairs.forEach((primerPair: PrimerPair) => {
            const color = validColors.pop();
            primerPair.color = color || '#ffffff';
          });
        });
        
        /* 
        Separate the primer pairs into single primers and append all into a single list.
        This way we can iterate over each primer in a pool if needed. This is useful for
        displaying all primers from a pool in a 'simple' way haha.
        */
        const primersByPool = new Map<string, Primer[]>();
        primerPairsGroupedByPoolId.forEach((primerPairs: PrimerPair[], key: string) => {
          // Extract the primers from the primer pairs
          const primers = primerPairs.reduce((acc: Primer[], element: PrimerPair) => {
            return acc.concat(element.forward_primer, element.reverse_primer);
          }, []);

          // Sort the primers by position
          primers.sort((a: Primer, b: Primer) => a.position - b.position);

          // Add the primers to the map
          primersByPool.set(key, primers);
        });
        
        // Recalculate start and end based on primers and apply a small buffer
        const start = Math.round(Math.min(...primerPairs.map((primerPair: PrimerPair) => primerPair.forward_primer.position))*0.99995);
        const end = Math.round(Math.max(...primerPairs.map((primerPair: PrimerPair) => primerPair.reverse_primer.position))*1.00005);
        // Create a new region object
        const region: Region = {
          name: regionInfo.name,
          start: start,
          end: end,
          sequence: regionInfo.sequence,
          primerPairsByPool: primerPairsGroupedByPoolId,
          primersByPool: primersByPool,
        };

        // Add the region to the list of regions
        results.set(key, region);
      });
    
    // Set the results
    this._results = results;
  }

  getResults(): Map<string, Region> {
    return this._results;
  }

  isErrored(): boolean {
    return this._error;
  }

  constructor() { 
    this._initResults();
  }
}
