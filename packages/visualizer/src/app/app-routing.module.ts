import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';
import { HomeComponent } from './home/home.component';
import { DetailedViewComponent } from './views/detailed-view/detailed-view.component';
import { SimplifiedViewComponent } from './views/simplified-view/simplified-view.component';

const routes: Routes = [
  { path: '', component: HomeComponent },
  { path: 'results/simple', component: SimplifiedViewComponent },
  { path: 'results/detail', component: DetailedViewComponent },
];

@NgModule({
  imports: [RouterModule.forRoot(routes)],
  exports: [RouterModule]
})
export class AppRoutingModule { }
