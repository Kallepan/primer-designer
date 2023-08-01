import { Injectable } from '@angular/core';
import { MatSnackBar } from '@angular/material/snack-bar';

const DURATION = 2000;

@Injectable({
  providedIn: 'root'
})
export class MessageService {

  constructor(
    private _snackBar: MatSnackBar,
  ) { }

  openSnackBar(message: string, action: string) {
    this._snackBar.open(message, action);
  }

  openTemporarySnackBar(message: string, action: string) {
    this._snackBar.open(message, action, { duration: DURATION });
  }
}
