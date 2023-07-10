import { Component, OnInit } from '@angular/core';
import { BehaviorSubject, tap } from 'rxjs';
import { ResultsService } from '../../services/results.service';
import { Region } from '../../types';
import { CONFIG } from '../../config';
import { MessageService } from 'src/app/services/message.service';

@Component({
  selector: 'app-detail-view',
  templateUrl: './detail-view.component.html',
  styleUrls: ['./detail-view.component.scss']
})
export class DetailViewComponent implements OnInit {
  minSequencePos = 1;
  maxSequencePos = 100;
  validRegions: string[] = [];
  currentSequence: string[] = [];
  presentPools: string[] = [];
  currentColors = new Map<string, string[]>();

  private _activeRegionData$ = new BehaviorSubject<Region | null>(null);
  activeRegionData$ = this._activeRegionData$.asObservable().pipe(
    tap(regionData => {
      // Check if the region data is available
      if (!regionData) return;

      // If regionData.sequence length is greater than the max sequence length display warning
      if (regionData.sequence.length > CONFIG.DETAILED_VIEW.SEQUENCE_MAX_LENGTH)
        this._messageService.openSnackBar(`Only the first ${CONFIG.DETAILED_VIEW.SEQUENCE_MAX_LENGTH} bases are displayed.`, 'OK');
      
      // Set the current sequence
      this.currentSequence = regionData.sequence.slice(
        0,
        CONFIG.DETAILED_VIEW.SEQUENCE_MAX_LENGTH,
      ).split('');

      // Set the min and max sequence positions
      this.minSequencePos = regionData.start;
      this.maxSequencePos = this.minSequencePos + this.currentSequence.length - 1;

      // Set the valid pools
      this.presentPools = Array.from(regionData.primersByPool.keys());

      // Update the color Encoding
      this._updateColorEncoding(regionData);
    }),
  );

  private _updateColorEncoding(regionData: Region): void {
    /* Updates the global color encoding for the region */

    // Reset the color encoding
    this.currentColors.clear();
    this.presentPools.forEach(pool_id => {
      // Create new color encoding for the pool
      const colors = regionData.sequence.split('').map((_, index) => {
        const color = this._getBgColor(index, pool_id);
        return color ? color : CONFIG.DETAILED_VIEW.DEFAULT_COLOR;
      });

      // Update the color encoding
      this.currentColors.set(pool_id, colors);
    });
  }

  // Get the color for the position in the pool
  getColor(position: number, pool_id: string): string {
    const colors = this.currentColors.get(pool_id);
    if (!colors) return CONFIG.DETAILED_VIEW.DEFAULT_COLOR;

    return colors[position];
  }
  // Get Tooltip for the position in the pool
  getTooltip(position: number, pool_id: string): string {
    const indexInRegion = position + this.minSequencePos;
    let tooltip = `Position: ${indexInRegion}`;

    const activeRegionData = this._activeRegionData$.value;
    if (!activeRegionData) return tooltip;

    // Check if position within an amplicon
    const amplicons = activeRegionData.primerPairsByPool.get(pool_id)?.filter(primerPair => {
      const primerStart = primerPair.forward_primer.position;
      const primerEnd = primerPair.reverse_primer.position + primerPair.reverse_primer.length;
      return indexInRegion >= primerStart && indexInRegion <= primerEnd;
    });
    // Check if position within a primer
    const primers = activeRegionData.primersByPool.get(pool_id)?.filter(primer => {
      const primerStart = primer.position;
      const primerEnd = primer.position + primer.length;
      return indexInRegion >= primerStart && indexInRegion <= primerEnd;
    });

    if (primers?.length) tooltip += `, Primer: ${primers[0].id}`;
    if (amplicons?.length) tooltip += `, Amplicon: ${amplicons[0].amplicon_name}`;

    return tooltip;
  }

  // Change the active region being displayed
  changeActiveRegion(event: any) {
    const regionName = event.value;
    if (!regionName) return;

    // Get the region from the results service
    const regionData = this._resultsService.getResults().get(regionName);
    if (!regionData) return;

    this._activeRegionData$.next(regionData);
  }

  // Color the background and letters based on primer adjecency
  private _getBgColor(position: number, pool_id: string): string {
    // Check if active region data is available
    const activeRegionData = this._activeRegionData$.value;
    if (!activeRegionData) return CONFIG.DETAILED_VIEW.DEFAULT_COLOR;

    // Get the index in the region
    const indexInRegion = position + this.minSequencePos;

    // Check if the position within an amplicon
    const amplicons = activeRegionData.primerPairsByPool.get(pool_id)?.filter(primerPair => {
      const ampliconStart = primerPair.forward_primer.position;
      const ampliconEnd = primerPair.reverse_primer.position + primerPair.reverse_primer.length;
      return indexInRegion >= ampliconStart && indexInRegion <= ampliconEnd;
    });
    // Check if position within a primer
    const primers = activeRegionData.primersByPool.get(pool_id)?.filter(primer => {
      const primerStart = primer.position;
      const primerEnd = primer.position + primer.length;
      return indexInRegion >= primerStart && indexInRegion <= primerEnd;
    });

    // Return the color prioritising primers over amplicons
    if (primers?.length) return CONFIG.DETAILED_VIEW.PRIMER_COLOR;
    if (amplicons?.length) return amplicons[0].color || CONFIG.DETAILED_VIEW.DEFAULT_COLOR;

    return CONFIG.DETAILED_VIEW.DEFAULT_COLOR;
  }

  ngOnInit(): void {
    // Initialize the valid regions
    const regionData = this._resultsService.getResults();
    this.validRegions = Array.from(regionData.keys());
  }

  constructor(
    private _resultsService: ResultsService,
    private _messageService: MessageService,
  ) { }
}
