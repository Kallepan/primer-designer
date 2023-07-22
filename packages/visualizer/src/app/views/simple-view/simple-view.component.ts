/*
  This component is responsible for displaying the simplified view of the results.
  Here all regions are 'plotted' as a line with coordinates. Above the region, the
  primer pairs are displayed as rectangles. The pools are represented as offsets of
  the primer pairs. The primer pairs are colored according to the pool they belong.
*/
import { Component, OnInit } from '@angular/core';
import { BehaviorSubject, tap, map } from 'rxjs';
import { ResultsService } from 'src/app/services/results.service';
import { Region } from 'src/app/types';
import { SimplifiedAmpliconData, SimplifiedPrimerData, SimplifiedRegionData } from './types';

@Component({
  selector: 'app-simple-view',
  templateUrl: './simple-view.component.html',
  styleUrls: ['./simple-view.component.scss']
})
export class SimpleViewComponent implements OnInit {
  loading = true;

  private _data$ = new BehaviorSubject<Map<string,Region> | null>(null);
  data$ = this._data$.asObservable().pipe(
    map(data => {
      // Check if data is available
      if (!data) return;

      // Format the data
      const formattedData: SimplifiedRegionData[] = [];

      // Iterate over each region and extract the position and primers
      data.forEach((regionData) => {
        // Fetch the start and end positions of the region
        const start = regionData.start;
        const end = regionData.end;

        // Iterate over each pool and extract the amplicons
        const formattedAmplicons: SimplifiedAmpliconData[] = [];
        regionData.primerPairsByPool.forEach((primerPairs, poolId) => {
          primerPairs.forEach(primerPair => {
            const ampliconName = primerPair.amplicon_name;
            const ampliconStart = primerPair.forward_primer.position + primerPair.forward_primer.length;
            const ampliconEnd = primerPair.reverse_primer.position;

            // Add the amplicon to the formatted data
            formattedAmplicons.push({
              name: ampliconName,
              pool: "Pool " + poolId,
              x1: ampliconStart,
              x2: ampliconEnd,
            });
          });
        });

        // Iterate over each pool and extract the primers
        const formattedPrimers: SimplifiedPrimerData[] = [];
        regionData.primersByPool.forEach((primers, poolId) => {
          primers.forEach(primer => {
            // Fetch start and end positions of the primer
            const primerStart = primer.position;
            const primerEnd = primer.position + primer.length;

            // Add the primer to the formatted data
            formattedPrimers.push({
              pool: "Pool " + poolId,
              id: primer.id,
              x1: primerStart,
              x2: primerEnd,
              sequence: primer.sequence,
            });
          });
        });

        const formattedRegionData: SimplifiedRegionData = {
          name: regionData.name,
          start: start,
          end: end,
          primers: formattedPrimers,
          amplicons: formattedAmplicons,
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
