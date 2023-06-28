import { ComponentFixture, TestBed } from '@angular/core/testing';

import { SimplifiedViewComponent } from './simplified-view.component';

describe('SimplifiedViewComponent', () => {
  let component: SimplifiedViewComponent;
  let fixture: ComponentFixture<SimplifiedViewComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [SimplifiedViewComponent]
    });
    fixture = TestBed.createComponent(SimplifiedViewComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
