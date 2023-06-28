import { Component, OnInit } from '@angular/core';
import { BehaviorSubject, tap } from 'rxjs';
import { ResultsService } from '../../services/results.service';
import { Region } from '../../types';
import { CONFIG } from '../../config';
import { MessageService } from 'src/app/services/message.service';

@Component({
  selector: 'app-detailed-view',
  templateUrl: './detailed-view.component.html',
  styleUrls: ['./detailed-view.component.scss']
})
export class DetailedViewComponent implements OnInit {
  minSequencePos = 0;
  maxSequencePos = 100;
  lowerPositionPointer = 0;
  upperPositionPointer = 100;
  validRegions: string[] = [];
  currentSequence: string[] = [];
  presentPools: string[] = [];
  positionColorEncodingPerPool = new Map<string, Map<number, string>>();

  private _activeRegionData$ = new BehaviorSubject<Region | null>(null);
  activeRegionData$ = this._activeRegionData$.asObservable().pipe(
    tap(regionData => {
      if (!regionData) return;

      // Update values when the active region changes
      this.minSequencePos = regionData.start;
      this.maxSequencePos = regionData.end;
      this.lowerPositionPointer = regionData.start;
      this.upperPositionPointer = regionData.start + CONFIG.DETAILED_VIEW.SEQUENCE_VISIBLE_LENGTH;

      const absolutePositionInString = 0;
      this.currentSequence = regionData.sequence.slice(
        absolutePositionInString,
        absolutePositionInString + CONFIG.DETAILED_VIEW.SEQUENCE_VISIBLE_LENGTH,
      ).split('');

      // Update the present pools
      this.presentPools = Array.from(regionData.primerPairsByPool.keys());
      this.presentPools.forEach(pool => {
        // Create a new Map for the pool
        const positionColorEncoding = new Map<number, string>();

        // Create batches of CONFIG.DETAILED_VIEW.BATCH_SIZE and color them by 'offloading' the work with Promises
        for (let i = 0; i < regionData.sequence.length; i += CONFIG.DETAILED_VIEW.BATCH_SIZE) {
          const batch = regionData.sequence.slice(i, i + CONFIG.DETAILED_VIEW.BATCH_SIZE);
          new Promise<null>(resolve => {
            batch.split('').forEach((letter, j) => {
              const bgColor = this.getBGColor(i + j + regionData.start, pool);
              if (!bgColor) {
                // This colors the background based on the DNA letter. Enabling this makes the plot difficult to read
                // positionColorEncoding.set(i + j + regionData.start, this.getDefaultBGColor(letter));
                return;
              }
              // Set the corresponsing color for the position
              positionColorEncoding.set(i + j + regionData.start, bgColor);
            });
            resolve(null);
          });
        }
        // Store the position color encoding for each pool globally
        this.positionColorEncodingPerPool.set(pool, positionColorEncoding);
      });
    }),
  );

  getTooltip(indexInRegion: number, pool_id: string): string {
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

  // Change the displayed sequence information
  onSliderInputChange() {
    // The string is indexed from 0 therefore the relative position is the absolute position - minSequencePos
    const relativeLowerPosition = this.lowerPositionPointer - this.minSequencePos;
    let relativeUpperPosition = this.upperPositionPointer - this.minSequencePos;

    // Check if the region data is available
    const regionData = this._activeRegionData$.value;
    if (!regionData) return;
    
    // Check if the difference between the lower and upper position is too large if so set the upper position to the config value
    if (relativeUpperPosition - relativeLowerPosition > CONFIG.DETAILED_VIEW.SEQUENCE_VISIBLE_LENGTH) {
      relativeUpperPosition = relativeLowerPosition + CONFIG.DETAILED_VIEW.SEQUENCE_VISIBLE_LENGTH;
      this._messageService.openSnackBar(`The maximum sequence length is ${CONFIG.DETAILED_VIEW.SEQUENCE_VISIBLE_LENGTH}`, 'OK');
    }
    
    // Update the sequence
    this.currentSequence = regionData.sequence.slice(
      relativeLowerPosition,
      relativeUpperPosition,
    ).split('');
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
  getBGColor(indexInRegion: number, pool_id: string): (string | undefined) {
    const activeRegionData = this._activeRegionData$.value;
    if (!activeRegionData) return undefined;

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
    if (amplicons?.length) return amplicons[0].color;

    return undefined;
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
