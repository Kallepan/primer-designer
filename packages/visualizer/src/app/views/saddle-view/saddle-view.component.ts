import { Component, OnInit } from '@angular/core';
import { BehaviorSubject, map, tap } from 'rxjs';
import { ResultsService } from 'src/app/services/results.service';
import { SaddlePlotData } from './types';

@Component({
  selector: 'app-saddle-view',
  templateUrl: './saddle-view.component.html',
  styleUrls: ['./saddle-view.component.scss']
})
export class SaddleViewComponent implements OnInit {
  loading = true;

  private _data$ = new BehaviorSubject<Map<string, Array<number>> | null>(null);
  data$ = this._data$.asObservable().pipe(
    map(data => {
      // Check if data is available
      if (!data) return;

      const formattedData: SaddlePlotData[] = [];
      data.forEach((losses, poolId) => {
        formattedData.push({
          name: poolId,
          losses: losses
        });
      });

      return formattedData;
    }),
    tap(data => {
      // Check if data is available
      if (!data) return;

      setTimeout(() => {
        this.loading = false;
      });
    })
  );

  ngOnInit(): void {
    // Fetch the data
    this._data$.next(this._resultsService.getLosses());
  }

  constructor(
    private _resultsService: ResultsService,
  ) { }
}
