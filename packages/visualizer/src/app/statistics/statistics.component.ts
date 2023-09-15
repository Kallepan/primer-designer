import { AfterViewInit, Component, OnInit, ViewChild } from '@angular/core';
import { MatTableDataSource } from '@angular/material/table';
import { ResultsService } from '../services/results.service';
import { FormattedAmplicon } from '../types';
import { MatPaginator } from '@angular/material/paginator';
import { MatSort } from '@angular/material/sort';

// Define the columns to display in the table
const COLUMNS = [
  //{key: 'region', label: 'Region'},
  { key: 'name', label: 'Amplicon' },
  { key: 'pool', label: 'Pool' },
  { key: 'discarded', label: 'Absent?' },
  { key: 'forwardPrimers', label: 'Total Forward Primers' },
  { key: 'discardedForwardPrimers', label: 'Discarded Forward Primers' },
  { key: 'reversePrimers', label: 'Total Reverse Primers' },
  { key: 'discardedReversePrimers', label: 'Discarded Reverse Primers' },
]

@Component({
  selector: 'app-statistics',
  templateUrl: './statistics.component.html',
  styleUrls: ['./statistics.component.scss']
})
export class StatisticsComponent implements OnInit, AfterViewInit {
  // Table columns
  displayedColumns = COLUMNS.map(c => c.key);
  columns = COLUMNS;

  // Pagination and sorting
  @ViewChild(MatPaginator) paginator: MatPaginator;
  @ViewChild(MatSort) sort: MatSort;

  totalDiscardedAmplicons: number = 0;
  totalPassedAmplicons: number = 0;
  totalAmplicons: number = this._resultsService.getAmplicons().length;

  // Get the amplicon data from the results service
  private _ampliconData = this._resultsService.getAmplicons().map(amplicon => {
    const formattedAmplicon: FormattedAmplicon = {
      name: amplicon.name,
      region: amplicon.region,
      pool: amplicon.pool,
      discarded: amplicon.discarded,
      forwardPrimers: `${amplicon.n_forward || "NA"}`,
      discardedForwardPrimers: `${amplicon.n_discarded_forward || "-"}`,
      reversePrimers: `${amplicon.n_reverse || "NA"}`,
      discardedReversePrimers: `${amplicon.n_discarded_reverse || "-"}`,
    };

    if (amplicon.discarded) {
      this.totalDiscardedAmplicons++;
    } else {
      this.totalPassedAmplicons++;
    }

    return formattedAmplicon;
  });

  dataSource = new MatTableDataSource(this._ampliconData);

  applyFilter(event: Event) {
    const value = (event.target as HTMLInputElement).value;
    this.dataSource.filter = value.trim().toLowerCase();
  }

  constructor(
    private _resultsService: ResultsService,
  ) { }

  ngOnInit(): void {
  }

  ngAfterViewInit(): void {
    this.dataSource.paginator = this.paginator;
    this.dataSource.sort = this.sort;
  }
}
