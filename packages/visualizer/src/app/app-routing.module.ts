import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';
import { HomeComponent } from './home/home.component';
import { DetailViewComponent } from './views/detail-view/detail-view.component';
import { SimpleViewComponent } from './views/simple-view/simple-view.component';
import { SaddleViewComponent } from './views/saddle-view/saddle-view.component';

const routes: Routes = [
  { path: '', component: HomeComponent },
  { path: 'results/simple', component: SimpleViewComponent },
  { path: 'results/detail', component: DetailViewComponent },
  { path: 'results/saddle', component: SaddleViewComponent },
  { path: '**', redirectTo: '' }
];

@NgModule({
  imports: [RouterModule.forRoot(routes)],
  exports: [RouterModule]
})
export class AppRoutingModule { }
