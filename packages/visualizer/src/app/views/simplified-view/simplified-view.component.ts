/*
  This component is responsible for displaying the simplified view of the results.
  Here all regions are 'plotted' as a line with coordinates. Above the region, the
  primer pairs are displayed as rectangles. The pools are represented as offsets of
  the primer pairs. The primer pairs are colored according to the pool they belong.
*/
import { Component, OnInit } from '@angular/core';
import { BehaviorSubject, tap, map } from 'rxjs';
import { CONFIG } from 'src/app/config';
import { ResultsService } from 'src/app/services/results.service';
import { Region } from 'src/app/types';

@Component({
  selector: 'app-simplified-view',
  templateUrl: './simplified-view.component.html',
  styleUrls: ['./simplified-view.component.scss']
})
export class SimplifiedViewComponent implements OnInit {
  regionNames: string[] = [];
  loading = true;
  
  private _data$ = new BehaviorSubject<Map<string,Region> | null>(null);
  data$ = this._data$.asObservable().pipe(
    tap(data => {
      if (!data) return;
      this.regionNames = Array.from(data.keys());
    }),
    map(data => {
      if (!data) return;
      
      type primerData = {
        pool: string;
        id: string;
        x1: number;
        x2: number;
        y: number;
      };

      const formattedData = new Map<string, primerData[]>();
      // Iterate over each region and calculate the relative positions for each primer
      data.forEach((regionData, regionName) => {
        // Fetch the start and end positions of the region
        const start = regionData.start;
        const end = regionData.end;
        const quotient = 1 / (end - start) * CONFIG.SIMPLIFIED_VIEW.SCALE_FACTOR;

        regionData.primersByPool.forEach((primers, poolId) => {
          const formattedPrimers: primerData[] = [];
          primers.forEach(primer => {
            // Fetch start and end positions of the primer
            const primerStart = primer.position;
            const primerEnd = primer.position + primer.length;

            // Relativize the start and end positions of the primer
            const newPrimerStart = (primerStart - start) * quotient;
            const newPrimerEnd = (primerEnd - start) * quotient;

            // Calculate y position of the primer by using the pool id and a constant
            const y = parseInt(poolId) * CONFIG.SIMPLIFIED_VIEW.PRIMER_OFFSET;

            // Add the primer to the formatted data
            formattedPrimers.push({
              pool: poolId,
              id: primer.id,
              x1: newPrimerStart,
              x2: newPrimerEnd,
              y: y,
            });
          });
          formattedData.set(regionName, formattedPrimers);
        });
      });
      setTimeout(() => this.loading = false, 0);
      console.log(formattedData);
      return formattedData;
    }),
  );

  ngOnInit(): void {
    // Initialize the region data and feed the BehaviorSubject
    const data = this._resultsService.getResults();
    this._data$.next(data);
  }

  constructor(
    private _resultsService: ResultsService,
  ) {  }
}
