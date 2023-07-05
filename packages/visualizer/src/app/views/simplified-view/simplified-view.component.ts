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
import { SimplifiedPrimerData, SimplifiedRegionData } from './types';

@Component({
  selector: 'app-simplified-view',
  templateUrl: './simplified-view.component.html',
  styleUrls: ['./simplified-view.component.scss']
})
export class SimplifiedViewComponent implements OnInit {
  loading = true;

  private _data$ = new BehaviorSubject<Map<string,Region> | null>(null);
  data$ = this._data$.asObservable().pipe(
    map(data => {
      // Check if data is available
      if (!data) return;

      // Format the data
      const formattedData: SimplifiedRegionData[] = [];

      // Iterate over each region and calculate the relative positions for each primer
      data.forEach((regionData) => {
        // Fetch the start and end positions of the region
        const start = regionData.start;
        const end = regionData.end;

        // quotient is used to relativize the start and end positions of the primer
        // It is calculated by dividing the scale factor by the length of the region
        //const quotient = 1 / (end - start) * CONFIG.SIMPLIFIED_VIEW.SCALE_FACTOR;

        // Iterate over each primer across all pool and append these 
        // to the formatted data
        
        const formattedPrimers: SimplifiedPrimerData[] = [];
        let idx = 1;
        regionData.primersByPool.forEach((primers, poolId) => {
          primers.forEach(primer => {
            // Fetch start and end positions of the primer
            const primerStart = primer.position;
            const primerEnd = primer.position + primer.length;

            // Relativize the start and end positions of the primer
            //const relativePrimerStart = (primerStart - start) * quotient;
            //const relativePrimerEnd = (primerEnd - start) * quotient;
            
            // Add the primer to the formatted data
            formattedPrimers.push({
              pool: "Pool " + poolId,
              id: primer.id,
              x1: primerStart,
              x2: primerEnd,
            });
          });

          // Increment the index after each pool
          idx += 1;
        });
        const formattedRegionData: SimplifiedRegionData = {
          name: regionData.name,
          start: start,
          end: end,
          primers: formattedPrimers,
        }

        // Add the formatted primers to the formatted data
        formattedData.push(formattedRegionData);
      });

      return formattedData;
    }),
    tap(data => {
      if (!data) return;
      
      // Set timeout to 0 to prevent ExpressionChangedAfterItHasBeenCheckedError
      setTimeout(() => {
        this.loading = false;
      });
    }),
  );

  ngOnInit(): void {
    this._data$.next(this._resultsService.getResults());
  }

  constructor(
    private _resultsService: ResultsService,
  ) {  }
}
